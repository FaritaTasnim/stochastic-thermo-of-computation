import imageio
from glob import glob

'''
Create separate file with desired png images and this file in it,
then run this file in IPython.
'''

files = glob('*.png')

images = []
for f in files:
    images.append(imageio.imread(f))

imageio.mimsave('out.gif', images, fps=5)
