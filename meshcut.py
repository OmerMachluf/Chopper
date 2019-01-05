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

INTERSECT_EDGE = 0
INTERSECT_VERTEX = 1
def distance(vertexA, vertexC, vertexD):
    distAC = (vertexA[0] - vertexC[0], vertexA[1] - vertexC[1], vertexA[2] - vertexC[2])
    distanceAC = math.sqrt(math.pow(distAC[0], 2) + math.pow(distAC[1], 2) + math.pow(distAC[2], 2))
    distAD = (vertexA[0] - vertexD[0], vertexA[1] - vertexD[1], vertexA[2] - vertexD[2])
    distanceAD = math.sqrt(math.pow(distAD[0], 2) + math.pow(distAD[1], 2) + math.pow(distAD[2], 2))
    return distanceAC < distanceAD

def split_model(mesh, plane):
    xBoundriesA = [0, 0]
    yBoundriesA = [0, 0]
    zBoundriesA = [0, 0]
    xBoundriesB = [0, 0]
    yBoundriesB = [0, 0]
    zBoundriesB = [0, 0]
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
                if mesh.verts[solovertice-1][0] > 0:
                    xBoundriesA = (0, mesh.verts[solovertice-1][0])
                if mesh.verts[solovertice-1][1] > 0:
                    yBoundriesA = (0, mesh.verts[solovertice-1][1])
                if mesh.verts[solovertice-1][2] > 0:
                    zBoundriesA = (0, mesh.verts[solovertice-1][2])

                if mesh.verts[solovertice-1][0] > 0:
                    xBoundriesA = (0, mesh.verts[solovertice-1][0])
                if mesh.verts[solovertice-1][1] > 0:
                    yBoundriesA = (0, mesh.verts[solovertice-1][1])
                if mesh.verts[solovertice-1][2] > 0:
                    zBoundriesA = (0, mesh.verts[solovertice-1][2])

                #calcs the boundries
                xBoundriesA[0] = min(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundriesA[1] = max(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundriesA[0] = min(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundriesA[1] = max(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundriesA[0] = min(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundriesA[1] = max(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartA.append(mesh.verts[solovertex-1])
                newverticesPartA.append(intersections[0][1])
                newverticesPartA.append(intersections[1][1])
                newvectorsA.append(mesh.normals[solovertex-1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsA.append(mesh.normals[solovertex-1])
                newvectorsA.append(mesh.normals[solovertex-1])

                #calcs the boundries
                xBoundriesB[0] = min(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundriesB[1] = max(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundriesB[0] = min(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundriesB[1] = max(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundriesB[0] = min(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundriesB[1] = max(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])

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
                #calcs the boundries
                xBoundriesB[0] = min(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundriesB[1] = max(mesh.verts[solovertex-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundriesB[0] = min(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundriesB[1] = max(mesh.verts[solovertex-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundriesB[0] = min(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundriesB[1] = max(mesh.verts[solovertex-1][2], intersections[0][1][2], intersections[1][1][2])

                newverticesPartB.append(mesh.verts[solovertex-1])
                newverticesPartB.append(intersections[0][1])
                newverticesPartB.append(intersections[1][1])
                newvectorsB.append(mesh.normals[solovertex-1])
            # Wrong normals - adding here simply to maintain track between numbers
                newvectorsB.append(mesh.normals[solovertex-1])
                newvectorsB.append(mesh.normals[solovertex-1])


                #calcs the boundries
                xBoundriesA[0] = min(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                xBoundriesA[1] = max(mesh.verts[groupvertex[0]-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                yBoundriesA[0] = min(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                yBoundriesA[1] = max(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                zBoundriesA[0] = min(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])
                zBoundriesA[1] = max(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])

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

                     #calcs the boundries
                    xBoundriesB[0] = min(mesh.verts[solovertex-1][0], )
                    xBoundriesB[1] = max(mesh.verts[solovertex-1][0], mesh.verts[groupvertex[1]-1][0], intersections[0][1][0], intersections[1][1][0])
                    yBoundriesB[0] = min(mesh.verts[solovertex-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                    yBoundriesB[1] = max(mesh.verts[groupvertex[0]-1][1], mesh.verts[groupvertex[1]-1][1], intersections[0][1][1], intersections[1][1][1])
                    zBoundriesB[0] = min(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])
                    zBoundriesB[1] = max(mesh.verts[groupvertex[0]-1][2], mesh.verts[groupvertex[1]-1][2], intersections[0][1][2], intersections[1][1][2])

                    newverticesPartB.append(mesh.verts[solovertex-1])  # Prevent duplicates possible vertices
                    newvectorsB.append(mesh.normals[solovertex-1])

                else:
                    keySolo = newverticesindexsB[solovertex]

                newverticesPartB.append(intersections[0][1]) # new vertices - how do we know normals?
                newverticesPartB.append(intersections[1][1])
                newvectorsB.append(mesh.normals[solovertex-1]) #wrong
                newvectorsB.append(mesh.normals[solovertex-1])

                if not oldvertices.has_key(groupvertex[0]):
                    keyGroup0 = len(newverticesPartA)
                    newverticesindexsA[groupvertex[0]] = keyGroup0
                    newverticesPartA.append(mesh.verts[groupvertex[0]-1])
                    newvectorsA.append(mesh.normals[groupvertex[0]-1])
                else:
                    keyGroup0 = newverticesindexsA[groupvertex[0]]
                
                if not oldvertices.has_key(groupvertex[1]):
                    keyGroup1 = len(newverticesPartA)
                    newverticesindexsA[groupvertex[1]] = keyGroup1
                    newverticesPartA.append(mesh.verts[groupvertex[1]-1])
                    newvectorsA.append(mesh.normals[groupvertex[1]-1]) 
                else:
                    keyGroup1 = newverticesindexsA[groupvertex[1]]
                
                if ((keyGroup0 == 141 and keyGroup1 == 143) or
                    (keyGroup1 == 141 and keyGroup0 == 143)):
                    print("break")
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
                    newverticesPartA.append(mesh.verts[solovertex-1])  # Prevent duplicates possible vertices
                    newvectorsA.append(mesh.normals[solovertex-1])
                else:
                    keySolo = newverticesindexsA[solovertex]

                newverticesPartA.append(intersections[0][1]) # new vertices - how do we know normals?
                newverticesPartA.append(intersections[1][1])

                newvectorsA.append(mesh.normals[solovertex-1]) #wrong
                newvectorsA.append(mesh.normals[solovertex-1])

                if not oldvertices.has_key(groupvertex[0]):
                    keyGroup0 = len(newverticesPartB)
                    newverticesindexsB[groupvertex[0]] = keyGroup0
                    newverticesPartB.append(mesh.verts[groupvertex[0]-1])
                    newvectorsB.append(mesh.normals[groupvertex[0]-1])
                else:
                    keyGroup0 = newverticesindexsB[groupvertex[0]]
                    
                if not oldvertices.has_key(groupvertex[1]):
                    keyGroup1 = len(newverticesPartB)
                    newverticesindexsB[groupvertex[1]] = keyGroup1
                    newverticesPartB.append(mesh.verts[groupvertex[1]-1]) 
                    newvectorsB.append(mesh.normals[groupvertex[1]-1])
                else:
                    keyGroup1 = newverticesindexsB[groupvertex[1]]
                
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
            newverticesPartB.append(mesh.verts[vertex])
            newvectorsB.append(mesh.normals[vertex])
        else:
            newverticesindexsA[vertex + 1] = len(newverticesPartA)
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
    modelB = TriangleMesh(newverticesPartB, newtrianglesPartB)
    modelB.normals = newvectorsB
    return (modelA, modelB)

def compute_vertex_side(plane, vertex, vertexAside, dist_tol=1e-8):
    """ True: Vertex on side B; False - Vertex on side A."""
    distA = point_to_plane_dist(vertexAside, plane)
    vertex = point_to_plane_dist(vertex, plane)
    return vertex * distA < 0


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
