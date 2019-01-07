import sys
import math
import pygame
import meshcut
import datetime
import numpy.linalg as la
from OpenGL.GL import *
from OpenGL.GLU import *
from pygame.constants import *
from pygame.locals import *
RED = (0,255,0)
WHITE = (255, 255, 255)
FILENAME = "Santa Claus2.obj"
b = 4
X_GOAL = 1
Y_GOAL = 2
Z_GOAL = 1

class TimedValue:

    def __init__(self):
        self._started_at = datetime.datetime.utcnow()

    def __call__(self):
        time_passed = datetime.datetime.utcnow() - self._started_at
        if time_passed.total_seconds() > 8:
            return False
        return True

class BSP:
    #BSP: [mesh1, mesh2, ...] - each mesh is a part of the model (all meshes together creates the whole model)
    def __init__(self, Object):
        self.meshs = []
        self.Object = Object

class OBJ:
    def __init__(self, filename, swapyz=False):
        # Loads a Wavefront OBJ file.
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []
        self.triangles = []
        self.mtl = None
        material = None
        for line in open(filename, "r"):
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue
            if values[0] == 'v':
                v = list(map(float, values[1:4]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.vertices.append(v)
            elif values[0] == 'vn':
                v = list(map(float, values[1:4]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.normals.append(v)
            elif values[0] == 'vt':
                self.texcoords.append(list(map(float, values[1:3])))
            elif values[0] in ('usemtl', 'usemat'):
                material = values[1]
            elif values[0] == 'mtllib':
                self.mtl = MTL(values[1])
            elif values[0] == 'f':
                face = []
                texcoords = []
                norms = []
                for v in values[1:4]:
                    w = v.split('/')
                    face.append(int(w[0]))
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                self.faces.append((face, norms, material))
                self.triangles.append(face)
        self.gl_list = glGenLists(1)
        glNewList(self.gl_list, GL_COMPILE)
        glEnable(GL_TEXTURE_2D)
        countNormals = 0
        countVertices = 0
        for face in self.faces:
            vertices, normals, material = face
            glBegin(GL_POLYGON)
            for i in range(len(vertices)):
                if normals[i] > 0:
                    countNormals = countNormals + 1
                    glNormal3fv(self.normals[normals[i] - 1])
                countVertices = countVertices + 1
                glVertex3fv(self.vertices[vertices[i] - 1])
            glEnd()
        glDisable(GL_TEXTURE_2D)
        glEndList()

def MTL(filename):
    contents = {}
    mtl = {}
    for line in open(filename, "r"):
        if line.startswith('#'): continue
        values = line.split()
        if not values: continue
        if values[0] == 'newmtl':
            mtl = contents[values[1]] = {}
        elif mtl is None:
            raise ValueError("mtl file doesn't start with newmtl stmt")
        elif values[0] == 'map_Kd':
            # load the texture referred to by this declaration
            mtl[values[0]] = values[1]
            surf = pygame.image.load(mtl['map_Kd'])
            image = pygame.image.tostring(surf, 'RGBA', 1)
            ix, iy = surf.get_rect().size
            texid = mtl['texture_Kd'] = glGenTextures(1)
            glBindTexture(GL_TEXTURE_2D, texid)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                            GL_LINEAR)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                            GL_LINEAR)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ix, iy, 0, GL_RGBA,
                         GL_UNSIGNED_BYTE, image)
        else:
            mtl[values[0]] = list(map(float, values[1:]))
    return contents

def EvalCuts(mesh, T2, P, newNormals, obj, beamCounter, counter = 0):
    resultSet = []
    orig = [0, 0, 0]
    for normal in newNormals:
        for i in range(2):
            if normal[0] > 0:
                orig[0] = orig[0] + 0.05
            elif normal[0] < 0:
                orig[0] = orig[0] - 0.05

            if normal[1] > 0:
                orig[1] = orig[1] + 0.05
            elif normal[1] < 0:
                orig[1] = orig[1] - 0.05

            if normal[2] > 0:
                orig[2] = orig[2] + 0.05
            elif normal[2] < 0:
                orig[2] = orig[2] - 0.05

            S = meshcut.calculateSPotential(P.normals, normal)
            plane = meshcut.Plane(tuple(orig), normal)
            counter = counter + 1
            modelA, modelB = meshcut.split_model(P, plane, S, beamCounter, counter)
            if (modelA != None and modelB != None):
                bspAfterSplit = BSP(mesh)
                meshsAfterSplit = []
                meshsAfterSplit.extend(T2)
                meshsAfterSplit.append(modelA)
                meshsAfterSplit.append(modelB)
                bspAfterSplit.meshs = meshsAfterSplit
                resultSet.append(bspAfterSplit)
                if (modelA.boundaries == False or modelB.boundaries == False):
                    print("false boundaries")
    return resultSet

def minSizes(BSPs, x_goal, y_goal, z_goal):
    minDiff = meshcut.calcDiffGoal(BSPs[0], x_goal, y_goal, z_goal)
    minBsp = BSPs[0]
    index = 0
    for i in range(len(BSPs)):
        diff = meshcut.calcDiffGoal(BSPs[i], x_goal, y_goal, z_goal)
        if (diff < minDiff):
            minDiff = diff
            minBsp = BSPs[i]
            index = i

    BSPs.pop(index)
    return minBsp


def HighestRanked(BSPs, mesh):
    maxRank = meshcut.f_util(BSPs[0], mesh)
    maxBSP = BSPs[0]
    index = 0
    for i in range(len(BSPs)):
        if (meshcut.f_util(BSPs[i], mesh) > maxRank):
            maxRank = meshcut.f_util(BSPs[i], mesh)
            maxBSP = BSPs[i]
            index = i

    BSPs.pop(index)
    return maxBSP

def beamSearch(mesh, b, x_goal, y_goal, z_goal, newNormals, obj):
    currentBSPs = []
    bsp = BSP(mesh)
    bsp.meshs.append(mesh)
    currentBSPs.append(bsp)
    beamCounter = 0
    while not (meshcut.allAtGoal(currentBSPs, x_goal, y_goal, z_goal)):
        newBSPs = []
        notAtGoalSet, currentBSPs = meshcut.notAtGoalSet(currentBSPs, x_goal, y_goal, z_goal)
        for bsp in notAtGoalSet:
            # T2 is meshs of the bsp T without P(the largest part in T)
            (P, T2) = meshcut.largestPart(bsp, x_goal, y_goal, z_goal) #p is a mesh
            beamCounter = beamCounter + 1
            newBSPs = newBSPs + EvalCuts(mesh, T2, P, newNormals, obj, beamCounter)
        while(len(newBSPs) > 0 and len(currentBSPs) < b):
            currentBSPs.append(minSizes(newBSPs, x_goal, y_goal, z_goal))
    if (len(currentBSPs) > 0):
        return HighestRanked(currentBSPs, mesh)
    else:
        return mesh

#main:
INTERSECT_EDGE = 0
INTERSECT_VERTEX = 1
pygame.init()
viewport = (800,600)
hx = viewport[0]/2
hy = viewport[1]/2
srf = pygame.display.set_mode(viewport, OPENGL | DOUBLEBUF)
glLightfv(GL_LIGHT0, GL_POSITION,  (-40, 200, 100, 0.0))
glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
glLightfv(GL_LIGHT0, GL_DIFFUSE, (0.5, 0.5, 0.5, 1.0))
glEnable(GL_LIGHT0)
glEnable(GL_LIGHTING)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_DEPTH_TEST)
glShadeModel(GL_SMOOTH)           # most obj files expect to be smooth-shaded
# LOAD OBJECT AFTER PYGAME INIT
obj = OBJ(FILENAME, swapyz=False)
orc = OBJ("octahedron2.obj", swapyz=False)
clock = pygame.time.Clock()
newNormals = []
hashDict = {}
for normal in orc.normals:
    sqr = math.pow(normal[0], 2) + math.pow(normal[1], 2) + math.pow(normal[2], 2)
    norm = math.sqrt(sqr)
    normalTuple = (normal[0] / norm, normal[1] / norm, normal[2] / norm)
    if (normalTuple in hashDict or
        (-normalTuple[0], -normalTuple[1], -normalTuple[2]) in hashDict):
        continue
    hashDict[normalTuple] = 1
    newNormals.append(normalTuple)
# TODO: need to add here: trying a few origins for each plane - need to add measurements of boundries of the model
# on a given vector direction so we can adjust the offset of the surface to be between the boundries
mesh = meshcut.TriangleMesh(tuple(obj.vertices), tuple(obj.triangles))
mesh.normals = obj.normals
#orig = [0, 0, 0]

choppedBsp = beamSearch(mesh, b, X_GOAL, Y_GOAL, Z_GOAL, newNormals, obj)
print(choppedBsp)

for mesh in choppedBsp.meshs:
    mesh.gl_list = glGenLists(1)
    glNewList(mesh.gl_list, GL_COMPILE)
    glEnable(GL_TEXTURE_2D)
    glFrontFace(GL_CCW)
    glBegin(GL_POLYGON)

    for triangle in mesh.tris:
        for i in range(len(triangle)):
            glNormal3fv(mesh.normals[triangle[i]])
            glVertex3fv(mesh.verts[triangle[i]])
    glEnd()
    glDisable(GL_TEXTURE_2D)
    glEndList()
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    width, height = viewport
    gluPerspective(90.0, width / float(height), 1, 100.0)
    glEnable(GL_DEPTH_TEST)
    glMatrixMode(GL_MODELVIEW)

    rx, ry = (0, 0)
    tx, ty = (0, 0)
    zpos = 5
    rotate = move = False
    value = TimedValue()
    while value():
        clock.tick(30)
        for e in pygame.event.get():
            if e.type == QUIT:
                sys.exit()
            elif e.type == KEYDOWN and e.key == K_ESCAPE:
                sys.exit()
            elif e.type == MOUSEBUTTONDOWN:
                if e.button == 4:
                    zpos = max(1, zpos - 1)
                elif e.button == 5:
                    zpos += 1
                elif e.button == 1:
                    rotate = True
                elif e.button == 3:
                    move = True
            elif e.type == MOUSEBUTTONUP:
                if e.button == 1:
                    rotate = False
                elif e.button == 3:
                    move = False
            elif e.type == MOUSEMOTION:
                i, j = e.rel
                if rotate:
                    rx += i
                    ry += j
                if move:
                    tx += i
                    ty -= j

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        # RENDER OBJECT
        glTranslate(tx / 20., ty / 20., - zpos)
        glRotate(ry, 1, 0, 0)
        glRotate(rx, 0, 1, 0)
        pygame.draw.ellipse(srf, RED, [300, 10, 50, 20])
        glCallList(mesh.gl_list)
        pygame.display.flip()