import math

def largestPart(BSP):
    largestMesh = None
    largestVolume = 0
    bspWithoutLargest = []
    for mesh in BSP.meshs:
        volume = VolumeOfMesh(mesh)
        if (volume > largestVolume):
            if (largestMesh != None):
                bspWithoutLargest.append(largestMesh)
            largestVolume = volume
            largestMesh = mesh
        else:
            bspWithoutLargest.append(mesh)
    return (largestMesh, bspWithoutLargest)

def allAtGoal(BSPs, x_goal, y_goal, z_goal):
    for BSP in BSPs:
        if (atGoal(BSP, x_goal, y_goal, z_goal) == False):
            return False
    return True

def notAtGoalSet(BSPs, x_goal, y_goal, z_goal):
    list = []
    for BSP in BSPs:
        if (atGoal(BSP, x_goal, y_goal, z_goal) == False):
            list.append(BSP)
    return list

def atGoal(BSP, x_goal, y_goal, z_goal):
    #if (len(BSP.meshs) == 0):
    #    return False
    for mesh in BSP.meshs:
        if (atMeshGoal(mesh, x_goal, y_goal, z_goal) == False):
            mesh.atGoal = False
            return False
        else:
            mesh.atGoal = True
    return True

def atMeshGoal(TriangleMesh, x_goal, y_goal, z_goal):
    if (TriangleMesh.boundaries == True):
        if (((TriangleMesh.xBoundaries[1] - TriangleMesh.xBoundaries[0]) > x_goal) or ((TriangleMesh.yBoundaries[1] - TriangleMesh.yBoundaries[0]) > y_goal) or (
                (TriangleMesh.zBoundaries[1] - TriangleMesh.zBoundaries[0]) > z_goal)):
            return False
        else:
            return True
    else:
        xBoundaries = [TriangleMesh.verts[0][0], TriangleMesh.verts[0][0]]
        yBoundaries = [TriangleMesh.verts[0][1], TriangleMesh.verts[0][1]]
        zBoundaries = [TriangleMesh.verts[0][2], TriangleMesh.verts[0][2]]
        for vertices in TriangleMesh.verts:
            if (vertices[0] < xBoundaries[0]):
                xBoundaries[0] = vertices[0]
            elif (vertices[0] > xBoundaries[1]):
                xBoundaries[1] = vertices[0]
            if (vertices[1] < yBoundaries[0]):
                yBoundaries[0] = vertices[1]
            elif (vertices[1] > yBoundaries[1]):
                yBoundaries[1] = vertices[1]
            if (vertices[2] < zBoundaries[0]):
                zBoundaries[0] = vertices[2]
            elif (vertices[2] > zBoundaries[1]):
                zBoundaries[1] = vertices[2]

        if (((xBoundaries[1]-xBoundaries[0]) > x_goal) or ((yBoundaries[1]-yBoundaries[0]) > y_goal) or ((zBoundaries[1]-zBoundaries[0]) > z_goal)):
            return False
        else:
            return True

def f_util(T, object):
    V = VolumeOfMesh(object)
    max = 0;
    for p in T.meshs:
        theta = calc_theta(p)
        Vp = VolumeOfMesh(p);
        m = 1-Vp/(theta*V)
        if (m > max):
            max = m

    return max

def calc_theta(p):
    return 10

#p1, p2, p3- vectors
def SignedVolumeOfTriangle(p1, p2, p3):
    #v321 = p3.X*p2.Y*p1.Z;
    #v231 = p2.X*p3.Y*p1.Z;
    #v312 = p3.X*p1.Y*p2.Z;
    #v132 = p1.X*p3.Y*p2.Z;
    #v213 = p2.X*p1.Y*p3.Z;
    #v123 = p1.X*p2.Y*p3.Z;
    v321 = p3[0] * p2[1] * p1[2];
    v231 = p2[0] * p3[1] * p1[2];
    v312 = p3[0] * p1[1] * p2[2];
    v132 = p1[0] * p3[1] * p2[2];
    v213 = p2[0] * p1[1] * p3[2];
    v123 = p1[0] * p2[1] * p3[2];

    return (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);


def VolumeOfMesh(mesh):
    vols = []
    for t in mesh.tris:
        t1 = t[0]
        t2 = t[1]
        t3 = t[2]
        p1 = mesh.verts[t1-1]
        p2 = mesh.verts[t2-1]
        p3 = mesh.verts[t3-1]
        #v1 = (p1.X-p2.X, p1.Y-p2.Y, p1.Z-p2.Z)
        #v2 = (p2.X-p3.X, p2.Y-p3.Y, p2.Z-p3.Z)
        #v3 = (p3.X-p1.X, p3.Y-p1.Y, p3.Z-p1.Z)
        v1 = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
        v2 = (p2[0] - p3[0], p2[1] - p3[1], p2[2] - p3[2])
        v3 = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])
        vols.append(SignedVolumeOfTriangle(v1, v2, v3))
    return math.fabs(sum(vols))