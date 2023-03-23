import numpy as np
import subprocess
import math
from decimal import Decimal
import re
import gmsh
import operator
import sys

inputfile = 'msh_to_nemesh.inp'

# ifbu = '0'
meshname_list = []
numlattice = 0
with open(inputfile, 'r') as readobj:
    lines = readobj.readlines()
    for line in lines:
        tmplist = line.split()
        if len(tmplist) > 0:
            firststr = tmplist[0]
            if firststr == 'LAT':
                numlattice += 1
                meshname_list.append(tmplist[1]) # 存储各个栅元的名称，对应产生XX.nemesh和XX.regionalias          

def changeform(oldlist):
    newlist = ['{:.5e}'.format(i) + '   ' for i in oldlist]
    return (newlist)


def listclean(oldlist):
    newlist = str(oldlist).replace('[', '').replace(']', '')
    newlist = newlist.replace("'", '').replace(',', '')
    return (newlist)


def round_off(tmplist):
    for i in range(len(tmplist)):
        if abs(tmplist[i]) < 1e-15:
            tmplist[i] = 0.0
        tmplist[i] = Decimal(tmplist[i]).quantize(Decimal("0.000000001"), rounding="ROUND_HALF_UP")


def strtoint(oldlist):
    testlist = [int(i) for i in oldlist]
    return (testlist)


def strtofloat(oldlist):
    testlist = [float(i) for i in oldlist]
    return (testlist)


def floattostr(oldlist):
    testlist = [str(i) for i in oldlist]
    return (testlist)

def iflisthasequal(list1, list2):
    tag = False
    for i in list1:
        if i in list2:
            tag = True
    
    return(tag)

def list_move_left(list, offset):
    for i in range(offset):
        list.insert(len(list), list[0])
        list.remove(list[0])
    # return templist
def list_move_right(list, offset):
    templist = list[:]
    for i in range(offset):
        templist.insert(0, templist.pop())
    return templist
    
# Element Type index in the gmsh
# 2:  3-node triangle.
# 3:  4-node quadrangle..
# 9:
# 6-node second order triangle (3 nodes associated with the vertices and
# 3 with the edges).
# 16:
# 8-node second order quadrangle (4 nodes associated with the vertices
# and 4 with the edges).
# Element type index in the VITAS_FE
# !  ElementType =   5   2-D   Triangular      Linear
# !  ElementType =   6   2-D   Triangular      Quadratic
# !  ElementType =  10   2-D   Quadrilateral   Linear
# !  ElementType =  11   2-D   Quadrilateral   Quadratic
# !-------------------------------------------------------------------------
# !  Element Vertex Layout in VITAS_FE:
# !   Y
# !   |  Order = 1       Order = 2
# !   1 3               5              
# !     |\              |\             
# !     |  \            |  \           
# !     |    \          6    4         
# !     |      \        |      \       
# !     |        \      |        \     
# !   0 1---------2     1----2----3    
# !     0         1     0         1  -> X  
# !  Element Vertex Layout in Gmsh:
# !   Y
# !   |  Order = 1       Order = 2
# !   1 3               3              
# !     |\              |\             
# !     |  \            |  \           
# !     |    \          6    5         
# !     |      \        |      \       
# !     |        \      |        \     
# !   0 1---------2     1----4----2    
# !     0         1     0         1  -> X 
# !  Element Vertex Layout in VITAS_FE:
# !   Y
# !   |    Order = 1               Order = 2
# !   1 4-------------3         7------6------5
# !     |             |         |             |
# !     |             |         |             |
# !     |             |         8             4
# !     |             |         |             |
# !     |             |         |             |
# !  -1 1-------------2         1------2------3
# !    -1             1        -1             1  -> X
# !  Element Vertex Layout in Gmsh:
# !   Y
# !   |    Order = 1               Order = 2
# !   1 4-------------3         4------7------3
# !     |             |         |             |
# !     |             |         |             |
# !     |             |         8             6
# !     |             |         |             |
# !     |             |         |             |
# !  -1 1-------------2         1------5------2
# !    -1             1        -1             1  -> X
# !-------------------------------------------------------------------------
# !-------------------------------------------------------------------------
# GMSH与VITAS_FE之间有限元类型编号的映射
def gmsh_to_vitas_fe(elementtype):
    typeindex_map = {
        '1' : '1',
        '2' : '5',
        '3' : '10',
        '4' : '15',
        '5' : '25',
        '6' : '20',
        '8' : '2',
        '9' : '6',
        '11' : '16',
        '16' : '11',
        '17' : '26',
        '18' : '21',
        '10' : '52'
    }
    return typeindex_map.get(elementtype, None)
# type = 0  X
# type = 1  Y  
# islargetosmall = True     from largest to smallest
# islargetosmall = False    from smallest to largest
# 使有限元节点按照逆时针排列
def reorder(globalnodetag, localnodecoord, type, islargetosmall):
    def takefloat(elem):
        return float(elem[type])
    # globalnodeindex = []
    # localnodecoord = []
    # newconnectivity = []
    # for boundnode in (boundnode_in_this_ele):
    #     for i in range(len(nodetag)):
    #         if nodetag[i] == boundnode:
    #             globalnodeindex.append(nodetag[i])
    #             localnodecoord.append(coordinate[i])
    oldnodeconnet = localnodecoord[:]
    returnlist = []
    sorted_localnodecoord = localnodecoord
    sorted_localnodecoord.sort(key=takefloat, reverse=islargetosmall)
    length = len(localnodecoord)
    for i in sorted_localnodecoord:
        for j in range(length):
            if oldnodeconnet[j] == i:
                returnlist.append(globalnodetag[j])
    return returnlist
# 使list向左移动offset
def list_remove_left(list, offset):
    templist = list[:]
    for i in range(offset):
        templist.insert(len(templist), templist[0])
        templist.remove(templist[0])
    return templist
# 使list向右移动offset
def list_remove_right(list, offset):
    templist = list[:]
    for i in range(offset):
        templist.insert(0, templist.pop())
    return templist
def listclean(oldlist):
    newlist = str(oldlist).replace('[', '').replace(']', '')
    newlist = newlist.replace("'", '').replace(',', '')
    return (newlist)
def ljust_list(oldlist, numspace):
    newlist = []
    for i in range(len(oldlist)):
        newlist.append(oldlist[i].ljust(numspace))
    return(newlist)

for meshindex in range(numlattice):
    gmshfile = meshname_list[meshindex]+'.msh'
    # 读取GMSH产生的.msh网格文件
    with open(gmshfile, 'r') as obj:
        lines = obj.readlines()
        for linenumber, line in enumerate(lines):
            temp = line.split()
            # 读取关键字所在的行号
            if '$PhysicalNames' in temp:
                PhysicalNamesline = linenumber
            if '$Entities' in temp:
                entityline = linenumber
            if '$EndEntities' in temp:
                endentityline = linenumber
                boundary_linenumber = linenumber - 1
            if '$Nodes' in temp:
                nodeline = linenumber
            if '$Elements' in temp:
                elementline = linenumber
            if '$EndElements' in temp:
                endeleline = linenumber
				
        # 读取physical group的信息
        numphysicalgroup = int(lines[PhysicalNamesline+1].split()[0])
        physicalgrouptag = [] # physical group tag
        physicalgroupname = [] # physical group name(即材料名称)
        matname = []
        for i in range(numphysicalgroup):
            physicalgrouptag.append(lines[PhysicalNamesline+2+i].split()[1])
            physicalgroupname.append(lines[PhysicalNamesline+2+i].split()[2])
            matname.append(lines[PhysicalNamesline+2+i].split()[2]) # 材料名称
			
        numpoint = int(lines[entityline + 1].split()[0]) # 几何点的数目
        numcurve = int(lines[entityline + 1].split()[1]) # 几何线的数目
        numsurface = int(lines[entityline + 1].split()[2]) # 几何面的数目
		
		# 读取每个几何面对应的physical group tag
        physicalgrouptag_of_surf = []
        for i in range(numsurface):
            temp = lines[entityline+2+numpoint+numcurve+i].split()
            physicalgrouptag_of_surf.append(temp[8]) 
			
        num_entities = endentityline - entityline         
        num_nodes = int(lines[nodeline + 1].split()[1]) # 有限元节点数目
        num_elements = int(lines[elementline + 1].split()[1]) # 有限元数目

        num_entities_block = int(lines[elementline + 1].split()[0])
        templine = nodeline + 2
        nodecoord = [] # 存储节点坐标
        coordline = 0
        boundarynode = []
        nodetag = []
        boundaryentity_dimtag = [] 

        for i in range(num_entities):
            if coordline < elementline - 1:
                temp = lines[templine].split()
                tempentitydim = int(temp[0])
                tempentitytag = int(temp[1])
                tempnumnode = int(temp[3])
                tempnodeline = templine + 1
                for node in range(tempnumnode):
                    tempnodetag = lines[tempnodeline].split()[0]
                    nodetag.append(tempnodetag) 
                    if tempentitydim == 0 and tempentitytag in boundnodetag:
                        boundaryentity_dimtag.append([tempentitydim, tempentitytag])
                        boundarynode.append(tempnodetag)
                    if tempentitydim == 1 and tempentitytag in boundlinetag:
                        boundaryentity_dimtag.append([tempentitydim, tempentitytag])
                        boundarynode.append(tempnodetag)
                    tempnodeline += 1
                coordline = templine + tempnumnode + 1
                for j in range(tempnumnode):
                    tempcoord = lines[coordline].split()
                    nodecoord.append(tempcoord[:2])
                    coordline += 1
                templine = coordline

        templine = elementline + 2
        node_in_element = []
        element_type_index = []
        mat_index = [] # 材料区编号
        elementtag = [] # 有限元编号
        ele_to_bureg = [] # 有限元编号向燃耗区编号的映射
        for i in range(num_entities_block):
            if templine < endeleline:
                temp = lines[templine].split()
                tempsurf = int(temp[1]) # 有限元对应的几何面的编号
                tempeletype = temp[2]
                tempnumelement = int(temp[3])
                node_in_ele_line = templine + 1 # 组成有限元的节点所在行的行号
                for j in range(tempnumelement):
                    element_type_index.append(tempeletype) # GMSH有限元类型编号
                    temp_phygrouptag = physicalgrouptag_of_surf[tempsurf-1]
                    temp_index = physicalgrouptag.index(temp_phygrouptag)
                    mat_index.append(temp_index+1) # 有限元对应的材料编号 1,2,3...

                    # if ifbu == '1':
                    #     ele_to_bureg.append(bureg_index[tempbrindex-1])
                    temp2 = lines[node_in_ele_line].split()
                    elementtag.append(temp2[0]) # 有限元编号
                    node_in_element.append(temp2[1:]) # 组成有限元的节点编号
                    node_in_ele_line += 1
                templine = node_in_ele_line
    eletypeindex_in_VITAS_FE = [] # 对应VITAS_FE中的有限元类型编号
    for i in range(num_elements):
        temptype = element_type_index[i]
        realtypeindex = gmsh_to_vitas_fe(temptype)
        eletypeindex_in_VITAS_FE.append(realtypeindex)
    numboundaryele = 0 # 边界上的有限元数目
    numboundary_surf = 0 # 位于边界上的有限元边的数目
    boundaryelement = [] # 存储边界上的有限元编号
    boundnode_in_ele = []
    elementorder = [] # 有限元阶数
    refsurf = [] # 参考面编号
    for i in range(num_elements):
        if eletypeindex_in_VITAS_FE[i] in ['6', '11']:
            temporder = 2 # 对应有限元阶数为2
        else:
            temporder = 1
        elementorder.append(temporder)
        numnode_in_boundary = 0 # 该有限元位于边界上的节点数目
        for j in node_in_element[i]:
            if j in boundarynode:
                boundnode_in_ele.append(j)
                numnode_in_boundary += 1
        if temporder == 1:
            # 当有限元阶数为1，若有限元多于两个节点在边界上，则该有限元属于边界有限元
            if numnode_in_boundary >= 2:
                numboundaryele += 1
                numboundary_surf += 1
                refsurf.append('2')
                boundaryelement.append(elementtag[i])
                if numnode_in_boundary > 2:
                    numboundary_surf += 1
                    boundaryelement.append(elementtag[i])
                    refsurf.append('3')
        else:
            # 当有限元阶数为2，若有限元多于3个节点在边界上，则该有限元属于边界有限元
            if numnode_in_boundary >= 3:
                numboundaryele += 1
                numboundary_surf += 1
                boundaryelement.append(elementtag[i])
                refsurf.append('2')
                if numnode_in_boundary > 3:
                    numboundary_surf += 1
                    boundaryelement.append(elementtag[i])  
                    refsurf.append('3')      
    # reorder the connectivity
    globaleletag = []
    boundaryeleindex = []
    for i in boundaryelement:
        for j in range(num_elements):
            if elementtag[j] == i:
                globaleletag.append(elementtag[j])
                boundaryeleindex.append(j)
    boundeleconnect = []
    for i in range(numboundary_surf):
        eleindex = boundaryeleindex[i]
        oldconnectivity = node_in_element[eleindex]
        boundnode_in_this_ele = []
        local_dimtag = []
        for j in range(len(boundarynode)):
            for k in range(len(oldconnectivity)):
                tempnodetag = oldconnectivity[k]
                if tempnodetag == boundarynode[j]:
                    boundnode_in_this_ele.append(tempnodetag)
                    local_dimtag.append(boundaryentity_dimtag[j])
        globalnodetag = []
        globalnodeindex = []
        localnodecoord = []
        newconnectivity = []
        for boundnode in (boundnode_in_this_ele):
            for ii in range(len(nodetag)):
                if nodetag[ii] == boundnode:
                    globalnodetag.append(nodetag[ii])
                    globalnodeindex.append(ii)
                    localnodecoord.append(nodecoord[ii])
                    kkkk = int(nodetag[ii])
                    if eletypeindex_in_VITAS_FE[eleindex] == '6': 
                        globalnodetag = globalnodetag[:3]
                    elif eletypeindex_in_VITAS_FE[eleindex] == '11':
                        globalnodetag = globalnodetag[:4]
        # def reorder(globalnodetag, localnodecoord, type, islargetosmall)
        # type = 0  X
        # type = 1  Y  
        # islargetosmall = True     from largest to smallest
        # islargetosmall = False    from smallest to largest    
    if iflisthasequal(bottomdimtag, local_dimtag) and iflisthasequal(leftdimtag, local_dimtag):
        sorted_connect = reorder(globalnodetag, localnodecoord, 1, True)  
    elif iflisthasequal(bottomdimtag, local_dimtag) and not(iflisthasequal(leftdimtag, local_dimtag)):
            sorted_connect = reorder(globalnodetag, localnodecoord, 0, False)         
    elif iflisthasequal(rightdimtag, local_dimtag) and iflisthasequal(bottomdimtag, local_dimtag):
        sorted_connect = reorder(globalnodetag, localnodecoord, 0, False)
    elif iflisthasequal(rightdimtag, local_dimtag) and not(iflisthasequal(bottomdimtag, local_dimtag)):
            sorted_connect = reorder(globalnodetag, localnodecoord, 1, False)
    elif iflisthasequal(topdimtag, local_dimtag) and iflisthasequal(rightdimtag, local_dimtag):
        sorted_connect = reorder(globalnodetag, localnodecoord, 1, False)
    elif iflisthasequal(topdimtag, local_dimtag) and not(iflisthasequal(rightdimtag, local_dimtag)):
            sorted_connect = reorder(globalnodetag, localnodecoord, 0, True)
    elif iflisthasequal(leftdimtag, local_dimtag) and iflisthasequal(topdimtag, local_dimtag):
        sorted_connect = reorder(globalnodetag, localnodecoord, 0, True)
    elif iflisthasequal(leftdimtag, local_dimtag) and not(iflisthasequal(topdimtag, local_dimtag)):
        sorted_connect = reorder(globalnodetag, localnodecoord, 1, True)
    # print(sorted_connect)
    tempconnectivity = oldconnectivity[:]
    for jj in range(1):
        for kk in range(len(oldconnectivity)):
            if oldconnectivity[kk] == sorted_connect[jj]:
                tempindex = kk
        if eletypeindex_in_VITAS_FE[eleindex] in ['5', '10']:
            # if tempconnectivity[tempindex - 1] == sorted_connect[jj + 1]:
            #     tempconnectivity.reverse()
            tempconnectivity.reverse()
            # print(tempconnectivity)
            while tempconnectivity[1] != sorted_connect[0]:
                list_move_left(tempconnectivity, 1)
            # print(tempconnectivity)
            newconnectivity = tempconnectivity
        # print(sorted_connect)
        tempconnetivity = oldconnectivity
        for jj in range(1):
            for kk in range(len(oldconnectivity)):
                if oldconnectivity[kk] == sorted_connect[jj]:
                    tempindex = kk
            if eletypeindex_in_VITAS_FE[eleindex] in ['5', '10']:
                if tempconnetivity[tempindex - 1] == sorted_connect[jj + 1]:
                    tempconnetivity.reverse()
                newconnectivity = list_remove_left(tempconnetivity, tempindex)
                newconnectivity = list_remove_right(newconnectivity, 1)
            else:
                if eletypeindex_in_VITAS_FE[eleindex] in ['6', '11']:
                    if eletypeindex_in_VITAS_FE[eleindex] == '6':
                        pointnodeconnect = tempconnetivity[:3]
                        linenodeconnect = tempconnetivity[3:]
                    if eletypeindex_in_VITAS_FE[eleindex] == '11':
                        pointnodeconnect = tempconnetivity[:4]
                        linenodeconnect = tempconnetivity[4:]
                    if pointnodeconnect[tempindex - 1] == sorted_connect[jj + 1]:
                        pointnodeconnect.reverse()
                        linenodeconnect.reverse()
                    newpointnodeconnect = list_remove_left(pointnodeconnect, tempindex)
                    newlinenodeconnect = list_remove_left(linenodeconnect, tempindex) 
                    newpointnodeconnect = list_remove_right(newpointnodeconnect, 1)
                    newlinenodeconnect = list_remove_right(newlinenodeconnect, 1)  
                    for iii in range(len(newpointnodeconnect)):
                        newconnectivity.append(newpointnodeconnect[iii])
                        newconnectivity.append(newlinenodeconnect[iii]) 
                    # newconnectivity = list_remove_right(newconnectivity, 2)
                    
        boundeleconnect.append(newconnectivity) # 存储所有有限元重排后的节点顺序
    # 输出VITAS_FE可以识别的有限元网格文件xx.nemesh
    nemesh_content = []
    nemesh_content.append('! ANL FINITE ELEMENT INPUT FILE DESCRIPTION  9 HEADER LINES ALWAYS')
    nemesh_content.append('! CARD TYPE 1:  (Input Style: 0-indexed 1-not indexed) (Debug Printing: 1-10)')
    nemesh_content.append('! CARD TYPE 2:  (# Elements) (# Nodes) (# Edit Regions) (# boundary element surfaces)') 
    nemesh_content.append('! CARD TYPE 3:  [Optional Index] (ElementType) (Material)  ! READ AS (ELEMENTTYPE(I),I=1,NUMELEMENTS)') 
    nemesh_content.append('! CARD TYPE 4:  [Optional Index] (Element Connectivity)    ! READ AS (CONNECTIVITY(J),J=1,ELEMENTVERTICES) per element') 
    nemesh_content.append('! CARD TYPE 5:  [Optional Index] (X) [Y] [Z]               ! READ AS (XYZ(I,J),J=1,NUMDIMENSIONS)          per mesh point') 
    nemesh_content.append('! CARD TYPE 6:  [Optional Index] (Element #) (Ref. Surf.) (bound. cond.) ! READ AS (BOUNDARYLIST(I,J),J=1,3)  per boundary element surface') 
    nemesh_content.append('! CARD TYPE 7:  [Optional Index] (Reaction rate) (# elements)  ! READ AS EDITREACTION(I),EDITREGELEMENTS(I)   per edit region') 
    nemesh_content.append('! CARD TYPE 8:  [Optional Index] (Edit Region Elements ! READ AS (EDITREGION(I,J),J=1,NUMELEMENTS) per edit region') 
    nemesh_content.append('0'.ljust(6) + '00'.ljust(6))
    nemesh_content.append(str(num_elements).ljust(6) + str(num_nodes).ljust(6) + '0'.ljust(6) + str(numboundary_surf).ljust(6)) 
    for i in range(num_elements):
        tempstr = str(i+1).ljust(6) + eletypeindex_in_VITAS_FE[i].ljust(6) + str(mat_index[i]).ljust(6)
        nemesh_content.append(tempstr) 
    for i in range(num_elements):
        tempconnectivity = []
        if elementtag[i] in boundaryelement:
            for j in range(numboundary_surf):
                if boundaryelement[j] == elementtag[i]:
                    tempconnectivity = boundeleconnect[j][:]
        else:
            tempconnectivity = node_in_element[i][:]
        # tempconnectivity = node_in_element[i][:]
        templist = ljust_list(tempconnectivity, 6)
        templist = listclean(templist)
        tempstr = str(i+1).ljust(6) + templist
        nemesh_content.append(tempstr)
        
    for i in range(num_nodes):
        templist = ljust_list(nodecoord[i], 30)
        templist = listclean(templist)
        tempstr = str(i+1).ljust(6) + templist
        nemesh_content.append(tempstr)
    for i in range(numboundary_surf):
        tempstr = str(i+1).ljust(6) + boundaryelement[i].ljust(6) + refsurf[i].ljust(6) + '0'.ljust(6)
        nemesh_content.append(tempstr)
    # for i in range(num_elements):
    #     tempstr = str(i+1).ljust(6) + elementtag[i].ljust(6) + str(ele_to_bureg[i]).ljust(6) 
    #     nemesh_content.append(tempstr)
    nemeshfile = meshname_list[meshindex]
    with open(nemeshfile, 'w') as write_obj:
        for i in range(len(nemesh_content)):
            write_obj.write(nemesh_content[i])
            write_obj.write('\n')

	# 产生对应的.regionalias文件
    regionaliasfile = meshname_list[meshindex] + '.regionalias'
    nummat = len(matname) # 栅元包含的材料数目
    regionalias = []
    regionalias.append(nemeshfile)
    for j in range(nummat):
        regionalias.append('ALIAS REGION_00000000' + str(j+1) + '  ' + 'R_' + matname[j])
    with open(regionaliasfile, 'w') as reobj:
        for k in range(len(regionalias)):
            reobj.write(regionalias[k] + '\n')                