import meshcut
import math
import sys
import math
import pygame
from OpenGL.arrays import *
from OpenGL_accelerate import *
import OpenGL_accelerate.arraydatatype
from OpenGL.GL import *
from OpenGL.GLU import *
from pygame.constants import *
from pygame.locals import *

#objective functions
def NumberOfPartsObjective(tree):
    return 0

def ConnectorFeasiblityObjective(tree):
    return 0

def StructObjective(tree):
    return 0

def FragileObjective(tree):
    return 0

def AestheticsObjective(tree):
    return 0

def Runfunction(function, newTree):
    return function


#general
def LargestPart(model, tree):
    return tree

def AddToBsp(tree, plane):
    return tree

def BeamSearch(tree, b):
    currentBSPS = {}

def IsDifferentEnough(newTree, resultTrees):
    return True

def HighestRanked(trees):
    return trees[0]

def AllAtGoal(trees):
    return True

def NotAtGoalSet(trees):
    return trees

def BeamSearch(model, b):
    currentBsps = []
    while not AllAtGoal(currentBsps):
        newBsps = []
        for tree in NotAtGoalSet(currentBsps):
            currentBsps.remove(tree)
            P = LargestPart(model, tree)
            newBsps.append(EvalCuts(tree, P))
        while count(currentBsps) < b:
            currentBsps.append(HighestRanked(newBsps))

    return HighestRanked(currentBsps)

def EvalCuts(tree, largestPart):
    objFunctionsResult = 0
    allFunctions = {}
    planes = {}
    NormalsUnion = {}
    currentTrees = {}
    resultTrees = {}
    for x in NormalsUnion:
        for function in allFunctions:
            for plane in planes: #intersecting P
                newTree = AddToBsp(tree, plane)
                currentTrees += newTree
                objFunctionsResult += Runfunction(function, newTree)

    for tree in currentTrees:
        if (IsDifferentEnough(newTree, resultTrees)):
            resultTrees += newTree

    return resultTrees
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

def LoadNormals(filename):
    model = OBJ(filename)
    newNormals = []
    hashDict = {}
    for normal in model.normals:
        normalTuple = (normal[0], normal[1], normal[2])
        if (normalTuple in hashDict):
            continue
        hashDict[normalTuple] = 1
        newNormals.append(normal)
    return newNormals

class OBJ:
    def __init__(self, filename, swapyz=False):
       # Loads a Wavefront OBJ file. 
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []

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
                for v in values[1:]:
                    w = v.split('/')
                    face.append(int(w[0]))
                    if len(w) >= 2 and len(w[1]) > 0:
                        texcoords.append(int(w[1]))
                    else:
                        texcoords.append(0)
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                self.faces.append((face, norms, texcoords, material)) 
        self  

normals = LoadNormals("C:\\Users\\omerm\\Downloads\\octahedron2.obj")
model = OBJ("C:\\3D Project\\Santa Claus.obj")
model.gl_list = glGenLists(1)
newVectors = []
vectors = mesh.get_attribute("vertex_normal")
mesh = meshcut.TriangleMesh(vectors, faces)
"""for vector in vectors:
    sqrsum = math.pow(vector.x, 2) + math.pow(vector.y, 2) + math.pow(vector.z, 2)
    sqrt = math.sqrt(sqrsum)
    newVector = (vector.x, vector.y, vector.z)
    if (sqrt != 1):
        newVector = (vector.x / sqrt, vector.y / sqrt, vector.z / sqrt)
    newVectors.append(newVector)
offset = (0, 0.5, 0)
for vector in newVectors:
    plane = meshcut.Plane((offset, vector))

textureSurface = pygame.image.load('pic.jpg')
textureData = pygame.image.tostring(textureSurface, "RGBA", 1)
hdc=windll.user32.GetDC(1)
print(hdc)
hglrc=wglCreateContext(hdc)
wglMakeCurrent(hdc,hglrc)
width = textureSurface.get_width()
height = textureSurface.get_height()
texture = glGenTextures(1)
glEnable(GL_TEXTURE_2D)
glBindTexture(GL_TEXTURE_1D, texture)
glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
    GL_UNSIGNED_BYTE, textureData)
model.gl_list = glGenLists(1) 
glNewList(model.gl_list, GL_COMPILE)
glEnable(GL_TEXTURE_2D)
glFrontFace(GL_CCW)
for face in self.faces:
    vertices, normals, texture_coords, material = face

    mtl = self.mtl[material]
    if 'texture_Kd' in mtl:
    # use diffuse texmap
        glBindTexture(GL_TEXTURE_2D, mtl['texture_Kd'])
    else:
    # just use diffuse colour
        glColor(*mtl['Kd'])

    glBegin(GL_POLYGON)
    for i in range(len(vertices)):
        if normals[i] > 0:
            glNormal3fv(self.normals[normals[i] - 1])
        if texture_coords[i] > 0:
            glTexCoord2fv(self.texcoords[texture_coords[i] - 1])
        glVertex3fv(self.vertices[vertices[i] - 1])
    glEnd()
glDisable(GL_TEXTURE_2D)
glEndList()
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

clock = pygame.time.Clock()

glMatrixMode(GL_PROJECTION)
glLoadIdentity()
width, height = viewport
gluPerspective(90.0, width/float(height), 1, 100.0)
glEnable(GL_DEPTH_TEST)
glMatrixMode(GL_MODELVIEW)

rx, ry = (0,0)
tx, ty = (0,0)
zpos = 5
rotate = move = False
while 1:
    clock.tick(30)
    for e in pygame.event.get():
        if e.type == QUIT:
            sys.exit()
        elif e.type == KEYDOWN and e.key == K_ESCAPE:
            sys.exit()
        elif e.type == MOUSEBUTTONDOWN:
            if e.button == 4: zpos = max(1, zpos-1)
            elif e.button == 5: zpos += 1
            elif e.button == 1: rotate = True
            elif e.button == 3: move = True
        elif e.type == MOUSEBUTTONUP:
            if e.button == 1: rotate = False
            elif e.button == 3: move = False
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
    glTranslate(tx/20., ty/20., - zpos)
    glRotate(ry, 1, 0, 0)
    glRotate(rx, 0, 1, 0)
    glCallList(model.gl_list)

    pygame.display.flip()"""