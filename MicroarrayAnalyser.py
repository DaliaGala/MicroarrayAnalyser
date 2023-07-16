#!/usr/bin/env pythonw

# Name: Microarray_analysis_program.py
# Created by: Dalia S. Gala, dalia.gala@merton.ox.ac.uk
# Date: 28 Oct 2019
# This program will accept a .jpg or .png file image containing multiple plates of DNA microarrays, and return a Comma Separated Variable (CSV) file listing the intensity of the signal strength for healthy (green) and unhealthy (red) test results.

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import signal
#from scipy.signal import find_peaks
import Pillow
from Pillow import Image
from Pillow import ImageFilter

# OPEN FILE AND DEFINE SUITABILITY
MicroarrayJpg = input("\033[1;37;41mEnter filename in current directory, or the entire filepath starting at home; your file HAS TO be a .jpg or a .png image file of a single symmetrical plate of a DNA microarray; if it is not, the program won\'t run correctly:\033[0m")

try:
	if len(MicroarrayJpg) == 0:
		print("Proper usage of this program requires the input of a filename of a .jpg or .png file")
		sys.exit() 
	elif ".jpg" in MicroarrayJpg or ".png" in MicroarrayJpg: 		#Input will require a .jpg file
		im = Image.open(MicroarrayJpg, 'r') 	#If input is a .jpg, program will run
	else:
		print("Invalid file: must be a .jpg file!") 
		sys.exit() 			#If input not jpg, error and exit
except FileNotFoundError:
	print("There is no such file - try again!")
	sys.exit() 				#If input is a .jpg but not exist, error and exit"""

#CONVERT TO BLACK AND WHITE
img = im
img = img.filter(ImageFilter.GaussianBlur(radius=2))
#img.show()  #DEBUGGING
img = img.convert('L') #makes it greyscale
#img.show()  #DEBUGGING
img_grey = np.asarray(img.getdata(),dtype=np.float64).reshape((img.size[1],img.size[0])) # get greyscale pixel values
#print(img_grey.shape)  #DEBUGGING
 
#GETTING AVERAGE VALUES FOR EACH COLUMN AND ROW; COMPARE UPON ROTATION TO GET MAXIMA
def GetAvgCols(imagefile):
	return(imagefile.mean(axis = 0)) #returns the average values of each column; axis 0 = means which list it reads
def GetAvgRows(imagefile):
	return(imagefile.mean(axis = 1)) # same as above, except for rows

GreyImgAvgRows = GetAvgRows(img_grey)
GreyImgAvgCols = GetAvgCols(img_grey)
#print(len(GreyImgAvgRows))  #DEBUGGING
#print(len(GreyImgAvgCols))  #DEBUGGING
#plt.plot(GreyImgAvgRows)  #DEBUGGING
#plt.show()  #DEBUGGING

# FIND PEAKS AND THEIR INDICES BLACK AND WHITE
# COLS
peaks_values_C = find_peaks(GreyImgAvgCols, distance = 10, prominence = 0.2, width = 2.5)

# ROWS
peaks_values_R = find_peaks(GreyImgAvgRows, distance = 10, prominence = 0.2, width = 2.5)

# SKIP THE 1ST ARRAY ELEMENT, ONLY TAKE 0TH
peaks_values_R = peaks_values_R[0]
peaks_values_C = peaks_values_C[0]
# print(peaks_values_R)  #DEBUGGING
# print(peaks_values_C)  #DEBUGGING

# plt.scatter(peaks_values_R, GreyImgAvgRows[peaks_values_R])  #DEBUGGING
# plt.show()  #DEBUGGING

# SORT INDICES BY SIZE
grid_rows_sorted = np.sort(peaks_values_R)
grid_cols_sorted = np.sort(peaks_values_C)
# sorts these by size to have the correct row and column order

# print(len(grid_rows_sorted)) #DEBUGGING
# print(len(grid_cols_sorted)) #DEBUGGING

# MAKE A GRID IMAGE WHICH WILL BE OVERLAID ON THE MICROARRAY TO SEE HOW WELL IT SHOWS WELLS
# where the lines cross - that is where the well centre should be - this section just to see the grid
width, height = im.size
plt.plot([])
plt.imshow(im)
for y in peaks_values_C:
	plt.axvline(width - y)
for x in peaks_values_R:
	plt.axhline(height - x)
cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_visible(False)
cur_axes.axes.get_yaxis().set_visible(False)
plt.savefig('GridforMicroarray.png')
# plt.show()  #DEBUGGING

# OPEN IMAGE IN COLOUR AND REFER THE PROGRAM TO THE COORDINATES OBTAINED ABOVE
# TO GET IT TO READ RGB VALUES
image = Image.open(MicroarrayJpg, 'r')
img = image.filter(ImageFilter.GaussianBlur(radius=2))
# img.show()  #DEBUGGING
rgb_image = image.convert('RGB')
r, g, b = rgb_image.getpixel((1, 1)) #converts the image to RGB and assigns the RGB values to r, g, b variables
# print(image)  #DEBUGGING

pixel_values = []
for coordR in range(len(peaks_values_R)):
	row_vec = []
	for coordC in range(len(peaks_values_C)):
		pixel = rgb_image.getpixel((int(peaks_values_C[coordC]), int(peaks_values_R[coordR])))
		# goes to the integer value of the pixel index in row and col and gets RGB values
		row_vec.append(pixel)
	pixel_values.append(row_vec)


pixel_values = np.array(pixel_values)
# print(pixel_values)  #DEBUGGING
# print(pixel_values.shape)  #DEBUGGING

# COMPARE R AND G VALUES AND APPEND THE DIFFERENCE TO AN ARRAY WHICH WILL BE PRINTED TO A CSV FILE
rg_values = []
for row in range(len(pixel_values)):
	rg_row = []
	for col in range(len(pixel_values[row])):
		if int(pixel_values[row][col][1]) - int(pixel_values[row][col][0]) > 0:
			rg_row.append(int(pixel_values[row][col][0]) - int(pixel_values[row][col][1]))
		elif int(pixel_values[row][col][0]) - int(pixel_values[row][col][1]) == 0:
			rg_row.append(0)
		else:
			rg_row.append(int(pixel_values[row][col][0]) - int(pixel_values[row][col][1]))
	rg_values.append(rg_row)
# "if" takes difference of G - R, if it is bigger than 0 (means more green), appends the value to row
# "elif" appends 0 if they are equal, so that should be a yellow well
# "else" appends the difference to row, and this time it will be negative because G - R is now <0
rg_values = np.array(rg_values)

#print(rg_values)  #DEBUGGING
#print(rg_values.shape)  #DEBUGGING

np.savetxt("well_values.csv", rg_values, delimiter=', ', newline='\n', fmt='%d')
#saves file as csv, gives new lines to each row, fmt expects integer values, %s = string, %f = float

