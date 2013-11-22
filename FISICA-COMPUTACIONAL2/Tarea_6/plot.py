# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
# <codecell>

#!/usr/bin/python
# -*- coding: latin-1 -*-
import math
import matplotlib.pyplot as plt
import numpy as np
import glob

"""
Programa que grafica archivos de texto con datos.
Autores: Oscar Fernando Garzon-Juan Diego Bermeo
Fecha: Martes 29 de Octubre de 2013.
"""
"""
Carga los datos desde un archivo "file" y los Grafica.
"""

# <codecell>

def Plot():
    """Plotea todos los datos sin importar el numero de archivos que entran como parametro"""

# <codecell>

"""Comenzamos cargando los datos."""

def Cargar_Datos():
    """Funcion que carga todos los archivos que estran como parametro para graficar"""
    


files = glob.glob('*.dat')
files
for i in range(len(files)):
    print files[-i-1]
files

for i in range(len(files)):
    archivo1=np.loadtxt(open(files[-i-1],"r"))
    
    for j in range (len(archivo1[:])):
       
        
        plt.scatter(archivo1[j][2],archivo1[j][1])
    plt.show()
    plt.close()
           



