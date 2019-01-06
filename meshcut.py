"""
Functions to slice a mesh. For now, computes planar cross-section
"""
import math
import numpy as np
import numpy.linalg as la
try:
    import scipy.spatial.distance as spdist
    USE_SCIPY = True
except ImportError:
    USE_SCIPY = False
import collections

# ---- Geometry datastructures


def make_edge(v1, v2):
    """
    We store edges as tuple where the vertex indices are sorted (so
    the edge going from v1 to v2 and v2 to v1 is the same)
    """
    tupsorted = tuple(sorted((v1, v2)))
    return tupsorted

def calculateSPotential(normals, normalT, THOLD = 0.95):
    S = []
    for i in range(len(normals)):
        dotProduct = normalT[0] * normals[i][0] + normalT[1] * normals[i][1] + normalT[2] * normals[i][2]
        if dotProduct > 0.95:
            S.append(i)
    return S

class TriangleMesh(object):
    def __init__(self, verts, tris):
        """
        Args:
            verts: The 3D vertex positions
            tris: A list of triplet containing vertex indices for each triangle
        """
        self.verts = np.array(verts)
        # For each edge, contains the list of triangles it b
        # ongs to
        # If the mesh is closed, each edge belongs to 2 triangles
        self.edges_to_tris = collections.defaultdict(lambda: [])
        # For each triangle, contains the edges it contains
        self.tris_to_edges = {}
        # For each vertex, the list of triangles it belongs to
        self.verts_to_tris = collections.defaultdict(lambda: [])
        self.boundaries = False
        self.tris = tris

        # Fill data structures
        for tid, f in enumerate(tris):
            if not isinstance(f, basestring):
                tri_edges = []
                for i in range(3):
                    v1 = f[i]
                    v2 = f[(i + 1) % 3]
                    e = make_edge(v1, v2)
                    if (self.edges_to_tris.has_key(e)):
                        val = self.edges_to_tris[e]
                        lst = list(val)
                        lst.append(tid)
                        self.edges_to_tris[e] = tuple(lst)
                    else:
                        self.edges_to_tris[e] = (tid,)
                    tri_edges.append(e)
                    if (self.verts_to_tris.has_key(f[i])):
                        val = self.verts_to_tris[f[i]]
                        lst = list(self.verts_to_tris[f[i]])
                        lst.append(tid)
                        self.verts_to_tris[f[i]] = tuple(lst)
                    else:
                        self.verts_to_tris[f[i]] = (tid,)
                self.tris_to_edges[tid] = tuple(tri_edges)
        
        # Sanity check : max 2 faces per edge
        for e, tris in self.edges_to_tris.items():
            if isinstance(tris, tuple):
                assert len(tris) <= 2

    def edges_for_triangle(self, tidx):
        """Returns the edges forming triangle with given index"""
        return self.tris_to_edges[tidx]

    def triangles_for_edge(self, edge):
        return self.edges_to_tris[edge]

    def triangles_for_vert(self, vidx):
        """Returns the triangles `vidx` belongs to"""
        return self.verts_to_tris[vidx]


class Plane(object):
    def __init__(self, orig, normal):
        self.orig = orig
        self.n = normal / la.norm(normal)

    def __str__(self):
        return 'plane(o=%s, n=%s)' % (self.orig, self.n)


def point_to_plane_dist(p, plane):
    return np.dot((p - plane.orig), plane.n)


def triangle_intersects_plane(mesh, tid, plane):
    """
    Returns true if the given triangle is cut by the plane. This will return
    false if a single vertex of the triangle lies on the plane
    """
    dists = [point_to_plane_dist(mesh.verts[vid], plane)
             for vid in mesh.tris[tid]]
    side = np.sign(dists)
    return not (side[0] == side[1] == side[2])


# ---- Planar cross-section
def triangle_contains_point(mesh, triangle, minX, maxX, minY, maxY, minZ, maxZ, point):
    if minX == maxX:
        if point[0] != minX:
            return False
        if point[1] > minY and point[1] < maxY and point[2] > minZ and point[2] < maxZ:
            return True
    if minY == maxY:
        if point[1] != minY:
            return False
        if point[0] > minX and point[0] < maxX and point[2] > minZ and point[2] < maxZ:
            return True
    if minZ == maxZ:
        if point[2] != minZ:
            return False
        if point[0] > minX and point[0] < maxX and point[1] > minY and point[1] < maxY:
            return True
    return False

INTERSECT_EDGE = 0
INTERSECT_VERTEX = 1
def distance(vertexA, vertexC, vertexD):
    distAC = (vertexA[0] - vertexC[0], vertexA[1] - vertexC[1], vertexA[2] - vertexC[2])
    distanceAC = math.sqrt(math.pow(distAC[0], 2) + math.pow(distAC[1], 2) + math.pow(distAC[2], 2))
    distAD = (vertexA[0] - vertexD[0], vertexA[1] - vertexD[1], vertexA[2] - vertexD[2])
    distanceAD = math.sqrt(math.pow(distAD[0], 2) + math.pow(distAD[1], 2) + math.pow(distAD[2], 2))
    return distanceAC < distanceAD

def split_model(mesh, plane, S, fragThreas = 0.1):
    S2 = []
    for vertex in S:
        dist = point_to_plane_dist(mesh.verts[vertex-1], plane)
        if dist < fragThreas:
            S2.append(vertex)
    xBoundariesA = [0, 0]
    yBoundariesA = [0, 0]
    zBoundariesA = [0, 0]
    xBoundariesB = [0, 0]
    yBoundariesB = [0, 0]
    zBoundariesB = [0, 0]
    allIntersections = cross_section_mesh(mesh, plane)
    if allIntersections == None or len(allIntersections) == 0:
        return (None, None)
    oldtriangles = {}
    oldvertices = {}
    newverticesindexsA = {}
    newverticesindexsB = {}
    newtrianglesPartA = []
    newtrianglesPartB = []
    newverticesPartA = []
    newverticesPartB = []
    newvectorsA = []
    newvectorsB = []
    indexA = 0
    indexB = 0
    # iterate through all intersections
    first = True
    for intersections in allIntersections:
        groupvertex = []
        # Starting with handling only 2 edges intersection. figuring out what do with vertexes later.
        if intersections[0][0] == INTERSECT_EDGE and intersections[1][0] == INTERSECT_EDGE:
            # Naive and stupid code to seperate sides -> 1 vertex vs 2
            edge0 = intersections[0][2]
            edge1 = intersections[1][2]
            if (edge0[0] == edge1[0] or edge0[0] == edge1[1]):
                solovertex = edge0[0]
            else:
                solovertex = edge0[1]
            if edge0[0] != solovertex:
                groupvertex.append(edge0[0])
            if edge0[1] != solovertex:
                groupvertex.append(edge0[1])
            if edge1[0] != solovertex:
                groupvertex.append(edge1[0])
            if edge1[1] != solovertex:
                groupvertex.append(edge1[1])
        else:
            continue

        tri = mesh.triangles_for_vert(solovertex) # get triangles
        trie0 = mesh.triangles_for_edge(edge0) # get the triangles for cutted edge 0
        trie1 = mesh.triangles_for_edge(edge1) # get the triangles for cutted ddge 1
        found = None
        for tri in trie0: # find the common triangle must be 1 only
            for tri1 in trie1:
                if tri == tri1:
                    found = tri
                    break
            if found != None:
                break

        oldtriangles[found] = 1
        tring = mesh.tris[found-1]
        ## Check if any point ofdffn S is fragile
        line = (intersections[0][1][0] - intersections[1][1][0], intersections[0][1][1] - intersections[1][1][1], intersections[0][1][2] - intersections[1][1][2])
        if (mesh.verts[tring[0]-1][0] == mesh.verts[tring[1]-1][0] and mesh.verts[tring[0]-1][0] == mesh.verts[tring[2]-1][0]):
            minX = maxX = mesh.verts[tring[0]-1][0]
            minY = min(mesh.verts[tring[0]-1][1], mesh.verts[tring[1]-1][1], mesh.verts[tring[2]-1][1])
            maxY = max(mesh.verts[tring[0]-1][1], mesh.verts[tring[1]-1][1], mesh.verts[tring[2]-1][1])
            minZ = min(mesh.verts[tring[0]-1][2], mesh.verts[tring[1]-1][2], mesh.verts[tring[2]-1][2])
            maxZ = max(mesh.verts[tring[0]-1][2], mesh.verts[tring[1]-1][2], mesh.verts[tring[2]-1][2])
        if (mesh.verts[tring[0]-1][1] == mesh.verts[tring[1]-1][1] and mesh.verts[tring[0]-1][10] == mesh.verts[tring[2]-1][1]):
            minY = maxY = mesh.verts[tring[0]-1][1]
            minX = min(mesh.verts[tring[0]-1][0], mesh.verts[tring[1]-1][0], mesh.verts[tring[2]-1][0])
            maxX = max(mesh.verts[tring[0]-1][0], mesh.verts[tring[1]-1][0], mesh.verts[tring[2]-1][0])
            minZ = min(mesh.verts[tring[0]-1][2], mesh.verts[tring[1]-1][2], mesh.verts[tring[2]-1][2])
            maxZ = max(mesh.verts[tring[0]-1][2], mesh.verts[tring[1]-1][2], mesh.verts[tring[2]-1][2])
        if (mesh.verts[tring[0]-1][2] == mesh.verts[tring[1]-1][2] and mesh.verts[tring[0]-1][2] == mesh.verts[tring[2]-1][2]):
            minZ = maxZ = mesh.verts[tring[0]-1][2]
            minX = min(mesh.verts[tring[0]-1][0], mesh.verts[tring[1]-1][0], mesh.verts[tring[2]-1][0])
            maxX = max(mesh.verts[tring[0]-1][0], mesh.verts[tring[1]-1][0], mesh.verts[tring[2]-1][0])
            minY = min(mesh.verts[tring[0]-1][1], mesh.verts[tring[1]-1][1], mesh.verts[tring[2]-1][1])
            maxY = max(mesh.verts[tring[0]-1][1], mesh.verts[tring[1]-1][1], mesh.verts[tring[2]-1][1])

        for s in S2:
            vert = mesh.verts[s-1]
            a = (vert[0] - intersections[0][1][0]) / line[0]
            b = (vert[1] - intersections[0][1][1]) / line[1]
            xA = intersections[0][1][0] + a * line[0]
            yA = intersections[0][1][1] + a * line[1]
            zA = intersections[0][1][2] + a * line[2]
            xB = intersections[0][1][0] + b * line[0]
            yB = intersections[0][1][1] + b * line[1]
            zB = intersections[0][1][2] + b * line[2]
            if vert[0] ==  xA and (vert[1] ==  yA or vert[2] == zA):
                pointonline = (xA, yA, zA)
            elif vert[0] ==  xB and (vert[1] ==  yB or vert[2] == zB):
                pointonline = (xB, yB, zB)
            if triangle_contains_point(mesh, tring, minX, maxX, minY, maxY, minZ, maxZ, pointonline):
                return (None, None)
        ### if we got here then no point is fragile.
        if (first):
            first = False
            xDist = pow(intersections[0][1][0] - intersections[1][1][0], 2)
            yDist = pow(intersections[0][1][1] - intersections[1][1][1], 2)
            zDist = pow(intersections[0][1][2] - intersections[1][1][2], 2)
            if (min(xDist, yDist, zDist) == xDist):
                compareX = 0
            elif (min(xDist, yDist, zDist) == yDist):
                compareX = 1
            elif (min(xDist, yDist, zDist) == zDist):
                compareX = 2   
            oldvertices[solovertex] = 1
            oldvertices[groupvertex[0]] = 1
            oldvertices[groupvertex[1]] = 1
            if mesh.verts[solovertex-1][compareX] > intersections[0][1][compareX]:
                #calcs the boundaries
                xBoundariesA[0] = min(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundariesA[1] = max(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundariesA[0] = min(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundariesA[1] = max(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundariesA[0] = min(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundariesA[1] = max(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartA.append(mesh.verts[solovertex-1])
                newverticesPartA.append(intersections[0][1])
                newverticesPartA.append(intersections[1][1])
                newvectorsA.append(mesh.normals[solovertex-1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsA.append(mesh.normals[solovertex-1])
                newvectorsA.append(mesh.normals[solovertex-1])

                #calcs the boundaries
                xBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartB.append(mesh.verts[groupvertex[0]-1])
                newverticesPartB.append(mesh.verts[groupvertex[1]-1])
                newvectorsB.append(mesh.normals[groupvertex[0]-1])
                newvectorsB.append(mesh.normals[groupvertex[1]-1])
                newverticesPartB.append(intersections[0][1])
                newverticesPartB.append(intersections[1][1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsB.append(mesh.normals[groupvertex[0]-1])
                newvectorsB.append(mesh.normals[groupvertex[1]-1])
            # TODO: 1. Add new normals for new vertices.
            #       2. Check the order - intersections[0] doesn't necessarily point to the intersection closest to groupvertex[0]
            #          If something is wrong here - expect to see self intersections. of triangles
                newtrianglesPartA.append((0, 1, 2))
                if distance(mesh.verts[groupvertex[0]-1], intersections[0][1], intersections[1][1]):
                    newtrianglesPartB.append((0, 1, 2))
                    newtrianglesPartB.append((1, 2, 3))
                else:
                    newtrianglesPartB.append((0, 1, 3))
                    newtrianglesPartB.append((1, 2, 3))
                newverticesindexsA[solovertex] = 0
                newverticesindexsB[groupvertex[0]] = 0
                newverticesindexsB[groupvertex[1]] = 1
                # Setting bounderies
                vertexA = solovertex
                vertexB = groupvertex[0]
            else:
                #calcs the boundaries
                xBoundariesB[0] = min(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundariesB[1] = max(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundariesB[0] = min(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundariesB[1] = max(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundariesB[0] = min(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundariesB[1] = max(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartB.append(mesh.verts[solovertex-1])
                newverticesPartB.append(intersections[0][1])
                newverticesPartB.append(intersections[1][1])
                newvectorsB.append(mesh.normals[solovertex-1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsB.append(mesh.normals[solovertex-1])
                newvectorsB.append(mesh.normals[solovertex-1])


                #calcs the boundaries
                xBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartA.append(mesh.verts[groupvertex[0]-1])
                newverticesPartA.append(mesh.verts[groupvertex[1]-1])
                newvectorsA.append(mesh.normals[groupvertex[0]-1])
                newvectorsA.append(mesh.normals[groupvertex[1]-1])
                newverticesPartA.append(intersections[0][1])
                newverticesPartA.append(intersections[1][1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsA.append(mesh.normals[groupvertex[0]-1])
                newvectorsA.append(mesh.normals[groupvertex[1]-1])
            # TODO: 1. Add new normals for new vertices.
            #       2. Check the order - intersections[0] doesn't necessarily point to the intersection closest to groupvertex[0]
            #          If something is wrong here - expect to see self intersections. of triangles
                newtrianglesPartB.append((0, 1, 2))
                if distance(mesh.verts[groupvertex[0]-1], intersections[0][1], intersections[1][1]):
                    newtrianglesPartA.append((0, 1, 2))
                    newtrianglesPartA.append((1, 2, 3))
                else:
                    newtrianglesPartA.append((0, 1, 3))
                    newtrianglesPartA.append((1, 2, 3))
                newverticesindexsB[solovertex] = 0
                newverticesindexsA[groupvertex[0]] = 0
                newverticesindexsA[groupvertex[1]] = 1
                # Setting bounderies
                vertexB = solovertex
                vertexA = groupvertex[0]
        else: # Seperate sides using original bounderies
            numofVerticesA = len(newverticesPartA)
            if (numofVerticesA == 141 or numofVerticesA == 143):
                print("break")
            if compute_vertex_side(plane, mesh.verts[solovertex-1], mesh.verts[vertexA-1]):
                if not oldvertices.has_key(solovertex):
                    keySolo = len(newverticesPartB)
                    newverticesindexsB[solovertex] = keySolo

                     #calcs the boundaries
                    xBoundariesB[0] = min(mesh.verts[solovertex-1][0], xBoundariesB[0])
                    xBoundariesB[1] = max(mesh.verts[solovertex-1][0], xBoundariesB[1])
                    yBoundariesB[0] = min(mesh.verts[solovertex-1][1], yBoundariesB[0])
                    yBoundariesB[1] = max(mesh.verts[solovertex-1][1], yBoundariesB[1])
                    zBoundariesB[0] = min(mesh.verts[solovertex-1][2], zBoundariesB[0])
                    zBoundariesB[1] = max(mesh.verts[solovertex-1][2], zBoundariesB[1])

                    newverticesPartB.append(mesh.verts[solovertex-1])  # Prevent duplicates possible vertices
                    newvectorsB.append(mesh.normals[solovertex-1])

                else:
                    keySolo = newverticesindexsB[solovertex]

                # calcs the boundaries
                xBoundariesB[0] = min(intersections[0][1][0], intersections[1][1][0], xBoundariesB[0])
                xBoundariesB[1] = max(intersections[0][1][0], intersections[1][1][0], xBoundariesB[1])
                yBoundariesB[0] = min(intersections[0][1][1], intersections[1][1][1], yBoundariesB[0])
                yBoundariesB[1] = max(intersections[0][1][1], intersections[1][1][1], yBoundariesB[1])
                zBoundariesB[0] = min(intersections[0][1][2], intersections[1][1][2], zBoundariesB[0])
                zBoundariesB[1] = max(intersections[0][1][2], intersections[1][1][2], zBoundariesB[1])

                newverticesPartB.append(intersections[0][1]) # new vertices - how do we know normals?
                newverticesPartB.append(intersections[1][1])
                newvectorsB.append(mesh.normals[solovertex-1]) #wrong
                newvectorsB.append(mesh.normals[solovertex-1])

                if not oldvertices.has_key(groupvertex[0]):
                    keyGroup0 = len(newverticesPartA)
                    newverticesindexsA[groupvertex[0]] = keyGroup0

                    # calcs the boundaries
                    xBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][0], xBoundariesA[0])
                    xBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][0], xBoundariesA[1])
                    yBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][1], yBoundariesA[0])
                    yBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][1], yBoundariesA[1])
                    zBoundariesA[0] = min(mesh.verts[groupvertex[0]-1][2], zBoundariesA[0])
                    zBoundariesA[1] = max(mesh.verts[groupvertex[0]-1][2], zBoundariesA[1])

                    newverticesPartA.append(mesh.verts[groupvertex[0]-1])
                    newvectorsA.append(mesh.normals[groupvertex[0]-1])
                else:
                    keyGroup0 = newverticesindexsA[groupvertex[0]]
                
                if not oldvertices.has_key(groupvertex[1]):
                    keyGroup1 = len(newverticesPartA)
                    newverticesindexsA[groupvertex[1]] = keyGroup1

                    # calcs the boundaries
                    xBoundariesA[0] = min(mesh.verts[groupvertex[1] - 1][0], xBoundariesA[0])
                    xBoundariesA[1] = max(mesh.verts[groupvertex[1] - 1][0], xBoundariesA[1])
                    yBoundariesA[0] = min(mesh.verts[groupvertex[1] - 1][1], yBoundariesA[0])
                    yBoundariesA[1] = max(mesh.verts[groupvertex[1] - 1][1], yBoundariesA[1])
                    zBoundariesA[0] = min(mesh.verts[groupvertex[1] - 1][2], zBoundariesA[0])
                    zBoundariesA[1] = max(mesh.verts[groupvertex[1] - 1][2], zBoundariesA[1])

                    newverticesPartA.append(mesh.verts[groupvertex[1]-1])
                    newvectorsA.append(mesh.normals[groupvertex[1]-1]) 
                else:
                    keyGroup1 = newverticesindexsA[groupvertex[1]]
                
                if ((keyGroup0 == 141 and keyGroup1 == 143) or
                    (keyGroup1 == 141 and keyGroup0 == 143)):
                    print("break")

                # calcs the boundaries
                xBoundariesA[0] = min(intersections[0][1][0], intersections[1][1][0], xBoundariesA[0])
                xBoundariesA[1] = max(intersections[0][1][0], intersections[1][1][0], xBoundariesA[1])
                yBoundariesA[0] = min(intersections[0][1][1], intersections[1][1][1], yBoundariesA[0])
                yBoundariesA[1] = max(intersections[0][1][1], intersections[1][1][1], yBoundariesA[1])
                zBoundariesA[0] = min(intersections[0][1][2], intersections[1][1][2], zBoundariesA[0])
                zBoundariesA[1] = max(intersections[0][1][2], intersections[1][1][2], zBoundariesA[1])
                newverticesPartA.append(intersections[0][1]) # new vertices- how do we know normals?
                newverticesPartA.append(intersections[1][1])

                newvectorsA.append(mesh.normals[groupvertex[0]-1]) #wrong
                newvectorsA.append(mesh.normals[groupvertex[1]-1])

                newtrianglesPartB.append((keySolo, len(newverticesPartB)-2, len(newverticesPartB)-1))
                if distance(mesh.verts[keyGroup0-1], intersections[0][1], intersections[1][1]):
                    newtrianglesPartA.append((keyGroup0, keyGroup1, len(newverticesPartA) - 2))
                    newtrianglesPartA.append((keyGroup1, len(newverticesPartA) - 2, len(newverticesPartA) - 1))
                else:
                    newtrianglesPartA.append((keyGroup0, keyGroup1, len(newverticesPartA) - 1))
                    newtrianglesPartA.append((keyGroup1, len(newverticesPartA) - 2, len(newverticesPartA) - 1))
            else: 
                if not oldvertices.has_key(solovertex):
                    keySolo = len(newverticesPartA)
                    newverticesindexsA[solovertex] = keySolo

                    # calcs the boundaries
                    xBoundariesA[0] = min(mesh.verts[solovertex-1][0], xBoundariesA[0])
                    xBoundariesA[1] = max(mesh.verts[solovertex-1][0], xBoundariesA[1])
                    yBoundariesA[0] = min(mesh.verts[solovertex-1][1], yBoundariesA[0])
                    yBoundariesA[1] = max(mesh.verts[solovertex-1][1], yBoundariesA[1])
                    zBoundariesA[0] = min(mesh.verts[solovertex-1][2], zBoundariesA[0])
                    zBoundariesA[1] = max(mesh.verts[solovertex-1][2], zBoundariesA[1])

                    newverticesPartA.append(mesh.verts[solovertex-1])  # Prevent duplicates possible vertices
                    newvectorsA.append(mesh.normals[solovertex-1])
                else:
                    keySolo = newverticesindexsA[solovertex]

                # calcs the boundaries
                xBoundariesA[0] = min(intersections[0][1][0], intersections[1][1][0], xBoundariesA[0])
                xBoundariesA[1] = max(intersections[0][1][0], intersections[1][1][0], xBoundariesA[1])
                yBoundariesA[0] = min(intersections[0][1][1], intersections[1][1][1], yBoundariesA[0])
                yBoundariesA[1] = max(intersections[0][1][1], intersections[1][1][1], yBoundariesA[1])
                zBoundariesA[0] = min(intersections[0][1][2], intersections[1][1][2], zBoundariesA[0])
                zBoundariesA[1] = max(intersections[0][1][2], intersections[1][1][2], zBoundariesA[1])

                newverticesPartA.append(intersections[0][1]) # new vertices - how do we know normals?
                newverticesPartA.append(intersections[1][1])

                newvectorsA.append(mesh.normals[solovertex-1]) #wrong
                newvectorsA.append(mesh.normals[solovertex-1])

                if not oldvertices.has_key(groupvertex[0]):
                    keyGroup0 = len(newverticesPartB)
                    newverticesindexsB[groupvertex[0]] = keyGroup0

                    # calcs the boundaries
                    xBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][0], xBoundariesB[0])
                    xBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][0], xBoundariesB[1])
                    yBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][1], yBoundariesB[0])
                    yBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][1], yBoundariesB[1])
                    zBoundariesB[0] = min(mesh.verts[groupvertex[0]-1][2], zBoundariesB[0])
                    zBoundariesB[1] = max(mesh.verts[groupvertex[0]-1][2], zBoundariesB[1])

                    newverticesPartB.append(mesh.verts[groupvertex[0]-1])
                    newvectorsB.append(mesh.normals[groupvertex[0]-1])
                else:
                    keyGroup0 = newverticesindexsB[groupvertex[0]]
                    
                if not oldvertices.has_key(groupvertex[1]):
                    keyGroup1 = len(newverticesPartB)
                    newverticesindexsB[groupvertex[1]] = keyGroup1

                    # calcs the boundaries
                    xBoundariesB[0] = min(mesh.verts[groupvertex[1]-1][0], xBoundariesB[0])
                    xBoundariesB[1] = max(mesh.verts[groupvertex[1]-1][0], xBoundariesB[1])
                    yBoundariesB[0] = min(mesh.verts[groupvertex[1]-1][1], yBoundariesB[0])
                    yBoundariesB[1] = max(mesh.verts[groupvertex[1]-1][1], yBoundariesB[1])
                    zBoundariesB[0] = min(mesh.verts[groupvertex[1]-1][2], zBoundariesB[0])
                    zBoundariesB[1] = max(mesh.verts[groupvertex[1]-1][2], zBoundariesB[1])

                    newverticesPartB.append(mesh.verts[groupvertex[1]-1])
                    newvectorsB.append(mesh.normals[groupvertex[1]-1])
                else:
                    keyGroup1 = newverticesindexsB[groupvertex[1]]

                # calcs the boundaries
                xBoundariesB[0] = min(intersections[0][1][0], intersections[1][1][0], xBoundariesB[0])
                xBoundariesB[1] = max(intersections[0][1][0], intersections[1][1][0], xBoundariesB[1])
                yBoundariesB[0] = min(intersections[0][1][1], intersections[1][1][1], yBoundariesB[0])
                yBoundariesB[1] = max(intersections[0][1][1], intersections[1][1][1], yBoundariesB[1])
                zBoundariesB[0] = min(intersections[0][1][2], intersections[1][1][2], zBoundariesB[0])
                zBoundariesB[1] = max(intersections[0][1][2], intersections[1][1][2], zBoundariesB[1])

                newverticesPartB.append(intersections[0][1]) # new vertices - how do we know normals?
                newverticesPartB.append(intersections[1][1])

                newvectorsB.append(mesh.normals[groupvertex[0]-1]) #wrong
                newvectorsB.append(mesh.normals[groupvertex[1]-1]) 

                newtrianglesPartA.append((keySolo, len(newverticesPartA)-2, len(newverticesPartA)-1))

                if distance(mesh.verts[keyGroup0-1], intersections[0][1], intersections[1][1]):
                    if (keyGroup0 == 91 or keyGroup1 == 91):
                        print("break")
                    newtrianglesPartB.append((keyGroup0, keyGroup1, len(newverticesPartB) - 2))
                    newtrianglesPartB.append((keyGroup1, len(newverticesPartB) - 2, len(newverticesPartB) - 1))
                else:
                    if (keyGroup0 == 91 or keyGroup1 == 91):
                        print("break")
                    newtrianglesPartB.append((keyGroup0, keyGroup1, len(newverticesPartB) - 1))
                    newtrianglesPartB.append((keyGroup1, len(newverticesPartB) - 2, len(newverticesPartB) - 1))
            oldvertices[solovertex] = 1
            oldvertices[groupvertex[0]] = 1
            oldvertices[groupvertex[1]] = 1
        # Now we are going to divide the common triangle to 3 triangles - 1 with 1 vertex and 2 intersection points(solo vertex)
        # 1 with 1 vertex and two intersection points(out of group vertex), and one with 2 vertexes and 1 intersection point(close to one not used)
        # Now create 3 lists: 1 for above planar cut model(add 1/2 triangles), 1 for below planar cut model(add rest of new triangles)
        # 1 for old triangles - so we'l know which ones to filter next.
        # iterate through all of the triangles and compose all untouched triangles - add to the right model
    newTrianglesForComparison = []
    newVerticesForComparison = []
    # compute additional vertices
    for i in range(len(mesh.verts)):
        if oldvertices.has_key(i+1):
            continue
        newVerticesForComparison.append(i)

    for vertex in newVerticesForComparison:
        if compute_vertex_side(plane, mesh.verts[vertex], mesh.verts[vertexA-1]):
            newverticesindexsB[vertex+1] = len(newverticesPartB)

            # calcs the boundaries
            xBoundariesB[0] = min(mesh.verts[vertex][0], xBoundariesB[0])
            xBoundariesB[1] = max(mesh.verts[vertex][0], xBoundariesB[1])
            yBoundariesB[0] = min(mesh.verts[vertex][1], yBoundariesB[0])
            yBoundariesB[1] = max(mesh.verts[vertex][1], yBoundariesB[1])
            zBoundariesB[0] = min(mesh.verts[vertex][2], zBoundariesB[0])
            zBoundariesB[1] = max(mesh.verts[vertex][2], zBoundariesB[1])

            newverticesPartB.append(mesh.verts[vertex])
            newvectorsB.append(mesh.normals[vertex])
        else:
            newverticesindexsA[vertex + 1] = len(newverticesPartA)

            # calcs the boundaries
            xBoundariesA[0] = min(mesh.verts[vertex][0], xBoundariesA[0])
            xBoundariesA[1] = max(mesh.verts[vertex][0], xBoundariesA[1])
            yBoundariesA[0] = min(mesh.verts[vertex][1], yBoundariesA[0])
            yBoundariesA[1] = max(mesh.verts[vertex][1], yBoundariesA[1])
            zBoundariesA[0] = min(mesh.verts[vertex][2], zBoundariesA[0])
            zBoundariesA[1] = max(mesh.verts[vertex][2], zBoundariesA[1])

            newverticesPartA.append(mesh.verts[vertex])
            newvectorsA.append(mesh.normals[vertex])
    # compute additional triangles
    for i in range(len(mesh.tris)):
        if (oldtriangles.has_key(i) or len(mesh.tris[i])>3):
            continue
        newTrianglesForComparison.append(mesh.tris[i])

    sideAtriangles, sideBtriangles = compute_triangle_sides(mesh, plane, newTrianglesForComparison, vertexA, vertexB, newverticesindexsA, newverticesindexsB)
    newtrianglesPartA = newtrianglesPartA + sideAtriangles
    newtrianglesPartB = newtrianglesPartB + sideBtriangles
    modelA = TriangleMesh(newverticesPartA, newtrianglesPartA)
    modelA.normals = newvectorsA
    modelA.boundaries = True
    modelA.xBoundaries = xBoundariesA
    modelA.yBoundaries = yBoundariesA
    modelA.zBoundaries = zBoundariesA

    modelB = TriangleMesh(newverticesPartB, newtrianglesPartB)
    modelB.normals = newvectorsB
    modelB.boundaries = True
    modelB.xBoundaries = xBoundariesB
    modelB.yBoundaries = yBoundariesB
    modelB.zBoundaries = zBoundariesB
    return (modelA, modelB)

def compute_vertex_side(plane, vertex, vertexAside, dist_tol=1e-8):
    """ True: Vertex on side B; False - Vertex on side A."""
    distA = point_to_plane_dist(vertexAside, plane)
    vertex = point_to_plane_dist(vertex, plane)
    return vertex * distA < 0

def largestPart(BSP):
    largestMesh = BSP.meshs[0]
    largestVolume = VolumeOfMesh(largestMesh)
    for mesh in BSP.meshs:
        volume = VolumeOfMesh(mesh)
        if (volume > largestVolume):
            largestVolume = volume
            largestMesh = mesh
    return largestMesh



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

def compute_triangle_sides(mesh, plane, triangles, vertexAside, vertexBside, dictA, dictB, dist_tol=1e-8):
    sideATriangles = []
    sideBTriangles = []
    newAVertices = []
    newBVertices = []
    distA = point_to_plane_dist(mesh.verts[vertexAside-1], plane)
    distB = point_to_plane_dist(mesh.verts[vertexBside-1], plane)
    for triangle in triangles:
        dists = []
        for vertice in triangle:
            dists.append(point_to_plane_dist(mesh.verts[vertice-1], plane))
        if (dists[0] * distA < 0 and
            dists[1] * distA < 0 and
            dists[2] * distA < 0):
            triangleN = (dictB[triangle[0]], dictB[triangle[1]], dictB[triangle[2]])
            sideBTriangles.append((triangleN))
        elif (dists[0] * distB < 0 and
            dists[1] * distB < 0 and
            dists[2] * distB < 0):
            triangleN = (dictA[triangle[0]], dictA[triangle[1]], dictA[triangle[2]])
            sideATriangles.append(triangleN)
    return (sideATriangles, sideBTriangles)

def compute_triangle_plane_intersections(mesh, tid, plane, dist_tol=1e-8):
    """
    Compute the intersection between a triangle and a plane

    Returns a list of intersections in the form
        (INTERSECT_EDGE, <intersection point>, <edge>) for edges intersection
        (INTERSECT_VERTEX, <intersection point>, <vertex index>) for vertices


    This return between 0 and 2 intersections :
    - 0 : the plane does not intersect the plane
    - 1 : one of the triangle's vertices lies on the plane (so it just
          "touches" the plane without really intersecting)
    - 2 : the plane slice the triangle in two parts (either vertex-edge,
          vertex-vertex or edge-edge)
    """
    # TODO: Use a distance cache
    dists = {}
    for vid in mesh.tris[tid]:
            if not isinstance(vid, basestring):
                dists[vid] = point_to_plane_dist(mesh.verts[vid-1], plane)
   # dists = {vid: point_to_plane_dist(mesh.verts[vid], plane)
    #         for vid in mesh.tris[tid]}
    # TODO: Use an edge intersection cache (we currently compute each edge
    # intersection twice : once for each tri)

    # This is to avoid registering the same vertex intersection twice
    # from two different edges
    vert_intersect = {vid: False for vid in dists.keys()}

    # Iterate through the edges, cutting the ones that intersect
    intersections = []
    
    edges = mesh.edges_for_triangle(tid)
    
    for e in edges:
        verts = []
        v1 = mesh.verts[e[0]-1]
        d1 = dists[e[0]]
        v2 = mesh.verts[e[1]-1]
        d2 = dists[e[1]]

        if np.fabs(d1) < dist_tol:
            # Avoid creating the vertex intersection twice
            if not vert_intersect[e[0]]:
                # point on plane
                intersections.append((INTERSECT_VERTEX, v1, e[0]))
                vert_intersect[e[0]] = True
        if np.fabs(d2) < dist_tol:
            if not vert_intersect[e[1]]:
                # point on plane
                intersections.append((INTERSECT_VERTEX, v2, e[1]))
                vert_intersect[e[1]] = True

        # If vertices are on opposite sides of the plane, we have an edge
        # intersection
        if d1 * d2 < 0:
            # Due to numerical accuracy, we could have both a vertex intersect
            # and an edge intersect on the same vertex, which is impossible
            if not vert_intersect[e[0]] and not vert_intersect[e[1]]:
                # intersection factor (between 0 and 1)
                # here is a nice drawing :
                # https://ravehgonen.files.wordpress.com/2013/02/slide8.png
                # keep in mind d1, d2 are *signed* distances (=> d1 - d2)
                s = d1 / (d1 - d2)
                vdir = v2 - v1
                ipos = v1 + vdir * s
                intersections.append((INTERSECT_EDGE, ipos, e))

    return intersections


def get_next_triangle(mesh, T, plane, intersection, dist_tol):
    """
    Returns the next triangle to visit given the intersection and
    the list of unvisited triangles (T)

    We look for a triangle that is cut by the plane (2 intersections) as
    opposed to one that only touch the plane (1 vertex intersection)
    """
    if intersection[0] == INTERSECT_EDGE:
        tris = mesh.triangles_for_edge(intersection[2])
    elif intersection[0] == INTERSECT_VERTEX:
        tris = mesh.triangles_for_vert(intersection[2])
    else:
        assert False, 'Invalid intersection[0] value : %d' % intersection[0]

    # Knowing where we come from is not enough. If an edge of the triangle
    # lies exactly on the plane, i.e. :
    #
    #   /t1\
    # -v1---v2-
    #   \t2/
    #
    # With v1, v2 being the vertices and t1, t2 being the triangles, then
    # if you just try to go to the next connected triangle that intersect,
    # you can visit v1 -> t1 -> v2 -> t2 -> v1 .
    # Therefore, we need to limit the new candidates to the set of unvisited
    # triangles and once we've visited a triangle and decided on a next one,
    # remove all the neighbors of the visited triangle so we don't come
    # back to it

    T = set(T)
    for tid in tris:
        if tid in T:
            intersections = compute_triangle_plane_intersections(
                    mesh, tid, plane, dist_tol)
            if len(intersections) == 2:
                T = T.difference(tris)
                return tid, intersections, T
    return None, [], T


def _walk_polyline(tid, intersect, T, mesh, plane, dist_tol):
    """
    Given an intersection, walk through the mesh triangles, computing
    intersection with the cut plane for each visited triangle and adding
    those intersection to a polyline.
    """
    T = set(T)
    p = []
    # Loop until we have explored all the triangles for the current
    # polyline
    while True:
        p.append(intersect[1])

        tid, intersections, T = get_next_triangle(mesh, T, plane,
                                                  intersect, dist_tol)
        if tid is None:
            break

        # get_next_triangle returns triangles that our plane actually
        # intersects (as opposed to touching only a single vertex),
        # hence the assert
        assert len(intersections) == 2

        # Of the two returned intersections, one should have the
        # intersection point equal to p[-1]
        if la.norm(intersections[0][1] - p[-1]) < dist_tol:
            intersect = intersections[1]
        else:
            assert la.norm(intersections[1][1] - p[-1]) < dist_tol, \
                '%s not close to %s' % (str(p[-1]), str(intersections))
            intersect = intersections[0]

    return p, T


def cross_section_mesh(mesh, plane, dist_tol=1e-8):
    """
    Args:
        mesh: A geom.TriangleMesh instance
        plane: The cut plane : geom.Plane instance
        dist_tol: If two points are closer than dist_tol, they are considered
                  the same
    """
    # Set of all triangles
    T = set(range(len(mesh.tris)))
    # List of all cross-section polylines
    P = []
    allIntersections = []
    while len(T) > 0:
        tid = T.pop()

        intersections = compute_triangle_plane_intersections(
                mesh, tid, plane, dist_tol)

        if len(intersections) == 2:
            allIntersections.append(intersections)
            #for intersection in intersections:
             #   p, T = _walk_polyline(tid, intersection, T, mesh, plane,
              #                        dist_tol)
               # if len(p) > 1:
                #    P.append(np.array(p))
    return allIntersections


def cross_section(verts, tris, plane_orig, plane_normal, **kwargs):
    """
    Compute the planar cross section of a mesh. This returns a set of
    polylines.

    Args:
        verts: Nx3 array of the vertices position
        faces: Nx3 array of the faces, containing vertex indices
        plane_orig: 3-vector indicating the plane origin
        plane_normal: 3-vector indicating the plane normal

    Returns:
        A list of Nx3 arrays, each representing a disconnected portion
        of the cross section as a polyline
    """
    mesh = TriangleMesh(verts, tris)
    plane = Plane(plane_orig, plane_normal)
    return cross_section_mesh(mesh, plane, **kwargs)


def pdist_squareformed_numpy(a):
    """
    Compute spatial distance using pure numpy
    (similar to scipy.spatial.distance.cdist())

    Thanks to Divakar Roy (@droyed) at stackoverflow.com

    Note this needs at least np.float64 precision!

    Returns: dist
    """
    a = np.array(a, dtype=np.float64)
    a_sumrows = np.einsum('ij,ij->i', a, a)
    dist = a_sumrows[:, None] + a_sumrows - 2 * np.dot(a, a.T)
    np.fill_diagonal(dist, 0)
    return dist


def merge_close_vertices(verts, faces, close_epsilon=1e-5):
    """
    Will merge vertices that are closer than close_epsilon.

    Warning, this has a O(n^2) memory usage because we compute the full
    vert-to-vert distance matrix. If you have a large mesh, might want
    to use some kind of spatial search structure like an octree or some fancy
    hashing scheme

    Returns: new_verts, new_faces
    """
    # Pairwise distance between verts
    if USE_SCIPY:
        D = spdist.cdist(verts, verts)
    else:
        D = np.sqrt(np.abs(pdist_squareformed_numpy(verts)))

    # Compute a mapping from old to new : for each input vert, store the index
    # of the new vert it will be merged into
    old2new = np.zeros(D.shape[0], dtype=np.int)
    # A mask indicating if a vertex has already been merged into another
    merged_verts = np.zeros(D.shape[0], dtype=np.bool)
    new_verts = []
    for i in range(D.shape[0]):
        if merged_verts[i]:
            continue
        else:
            # The vertices that will be merged into this one
            merged = np.flatnonzero(D[i, :] < close_epsilon)
            old2new[merged] = len(new_verts)
            new_verts.append(verts[i])
            merged_verts[merged] = True

    new_verts = np.array(new_verts)

    # Recompute face indices to index in new_verts
    new_faces = np.zeros((len(faces), 3), dtype=np.int)
    for i, f in enumerate(faces):
        new_faces[i] = (old2new[f[0]], old2new[f[1]], old2new[f[2]])

    # again, plot with utils.trimesh3d(new_verts, new_faces)
    return new_verts, new_faces

def f_util(T, object):
    V = VolumeOfMesh(object)
    max = 0;
    for p in T:
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

