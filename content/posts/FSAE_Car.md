---
title: "FSAE_Car"
date: 2022-02-08T15:12:25+05:30
draft: true

categories: [Converge 3.0]
tags: [CFD, Aerodynamics]
toc: false
author: "Satya Sasank"
---
AIM:

The aim of this project is to simulate the flow around an FSAE car in converge.

OBJECTIVES OF THE PROJECT:

Clean up the geometry     
Flag the parts as separate boundaries
Create a virtual wind tunnel
Run the simulation for two different types of Races

Race-1 

Racing track conditions - Lot of turns
Average lap speeds - 45kmph
70% turns at 45 degrees
20% turns at 80 degrees
10% turns at 20 degrees

Race-2

Racing track conditions – Straights
Average lap speed - 75kmph

PROCEDURE:

The Formula SAE is a student design competition organized by SAE international. Where students conceive, design, fabricate, and compete with small formula-style racing cars.

The car for our project is this [one](https://drive.google.com/drive/folders/1RHH6qC3BRMcH_rkYizZMvQxWGwteg06b?usp=sharing).

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (176)_1620966605.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (177) - Copy_1620968213.png) 

We import the STL file into Converge.

Once the geometry appears, we need to check the geometry units. We can do this by going to the Left-side geometry panel > Options > check the boundary box.

Converge uses the MKS system, while the geometry is created in IPS, so we need to transform the body to an appropriate scale to get meaningful results.

We can scale it by geometry panel > Transform > Scale > give a scaling factor of 0.0254. Before we create the wind tunnel, we need to know the dimensions of the wind tunnel and also we need to flag the boundaries of the car.

We use the boundary fence to construct boundary fences around parts that are to be flagged as separate boundaries.

And we flag the boundaries as

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (181)_1620976282.png)

Now to create the wind tunnel. We know that the wind tunnel must have at least 5 times the characteristic length in front of the body and at least 10 times the length behind the body. A height of at least twice the characteristic length above and beside the body.

With this information, we can figure out the wind tunnel dimensions based on the characteristic length which is 2.5 m (approx)

The calculated wind tunnel dimensions are as follows

X = 32.5 m

Y = 5 m

z = 2.5 m

To create the wind tunnel, We first create the vertices. For the vertices, we take the bottom of the tires so that the bottom of the tires should align with the bottom of the wind tunnel.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/InkedScreenshot (182)_LI_1620982074.jpg)

In the above picture, we deleted the bottom portion of the tires using the angle method.

All the points on the red circle are lying on the same plane. This plane which is parallel to the Z plane needs to be aligned with the bottom of the wind tunnel.

We first create the vertices for the wind tunnel using the wind tunnel dimensions and the location of the tire bottom and then patch these vertices using the ordered vertex method to form a wind tunnel.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (183)_1620983200.png)

Once the wind tunnel is created we flag the inlet and outlet boundaries. Now once we check for errors we should see all green in the diagnosis toolbar.

For different race conditions, we need to tilt the model in different angles as follows

20 deg
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (193)_1620999424.png)

45 deg
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (194)_1620999441.png)

80 deg
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (195)_1620999456.png)

We have to be careful while rotating the models because the triangles below the models tend to overlap. When this happens we need to delete the triangles that are overlapped and reconstruct those triangles. After reconstructing we need to check the diagnosis toolbar again for errors. Once all the checks are green we can follow the below procedure for different race tracks.

We then move onto case setup and then give the inputs as follows

Application: Time based

Gas simulations, Global transport properties remain the default values. N2 and O2 are the species to be added.

In the Run parameters, we select a transient solver. This is because the converge's steady-state solver gives inaccurate and gives results that are physically not possible.

the settings will be as follows
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (184)_1620991121.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (185)_1620991178.png)


All the boundaries will be given the Law of wall condition.

except for the wind tunnel. The wind tunnel inlet will have the same velocity as the region
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (186)_1620991394.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (187)_1620991473.png)


The above conditions are for race 2 for race 1 the velocity will be changed to 12.5m/s in both inlet and region 0.

The bottom will have the Law of the wall Boundary condition. All the boundaries are assigned to the same region 0

The turbulence model is STANDARD K Omega 2006. And since we cant achieve a cell size suitable for a y+ of 30-100 we use automatic wall function for wall treatment to obtain reasonable results at the walls.

In the Grid, control select Fixed embedding, Grid scaling options, and press ok. Two options appear in grid control.

In the Grid Control > Base Grid > set the x, y, z sizes to 0.3m in all the directions and click ok. with the fixed embedding enabled this is the closest we can get without crossing the 500000 cell limit.

For the fixed embedding we insert two types of embeddings one for the boundaries and another box type
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (188)_1620993371.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (189)_1620993388.png)

We give the same embedding for all the boundaries.

The grid scaling will be given based on a file and the inputs are as follows
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (190)_1620993816.png)


We then set the output file settings as follows. we can choose to select any additional variable that we need in post variable selection such as vorticity and set the output files as follows
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (192)_1620994083.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (191)_1620994067.png)


Once all the setup is done we do the final validation. Once every tab is green we export the files using the export input files option.  

We then run the simulation using the CYGWIN terminal. In the terminal, we navigate to the folder that contains the input files enter the following command.

"mpiexec.exe -n 4 converge.exe restricted </dev/null> logfile.txt &"

mpiexec.exe is the message passing interface executive file.
-n 6 represents the number of cores to be used for running the simulation.
converge.exe is the converge file that runs the simulation.
restricted limits the total number of cells to 5,00,000
</dev/null> logfile.txt is to create a log file that can be used to track the simulation.
& To keep the Cygwin terminal active for other commands

"tail -f logfile.txt"

This command tracks the log file and displays the content of the log file in the terminal.

The simulation might take several hours based on parameters like Hardware, grid size, time for which the simulation needs to be run. Once the calculations are finished exit the Cygwin terminal.

RESULTS AND POST-PROCESSING:

In the input file directory, a new folder will be created with the name "output". The post output files are stored here. We navigate to the post-processing module in Converge Studio and select this output and convert them into "Paraview-inline binary format" or "Ensight" format. 

Once this is done we open the results in paraview. We reduce the opacity of the side walls to see FSAE car inside the wind tunnel.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (207)_1621017302.png)

We can see that each boundary has been assigned with different colors. We then take a slice of this with Y-axis as the normal. This would cut the FSAE car in symmetric and we can observe the mesh using the surface-with-edges option. The below pictures are captured using the clip option since paraview was crashing upon slicing the geometry.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (204)_1621017323.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (205)_1621017336.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (206)_1621017358.png)
 

Instead of solid color, we can set it to velocity and pressure to see contours like these

Velocity
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (197)_1621017426.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (198)_1621017441.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (199)_1621017460.png)

Pressure
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (200)_1621017493.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (201)_1621017513.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (202)_1621017534.png)

Y plus
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Screenshot (203)_1621017560.png)

If we use the converge's line plotting module to plot the lift and drag for all the parts it would show up like this

Lift Force
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/LiftForce_1621014906.png)

Drag Force
![]https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Dragforce_1621014942.png)
         

If we export the same data into excel and look at it using python

(Code requires NumPy, pandas, and matplotlib libraries)               
```python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style

#Loading the workbook and reading the column headings excluding the first column
df = pd.read_excel('Lift.xlsx', index_col=0)
df_1 = pd.read_excel('Drag.xlsx', index_col=0)


#storing the column headings in a variable
components = df.columns
drag_components = df_1.columns


#declaring empty array for all the means
means = []
means_1 = []

for i in range (len(components)):
    #declaring empty array x  
    x = [] 
    
    #for i = 1; length (components) take the values of all rows in components this appends till length(components)
    x.append(df[components[i]])
    #Averaging the values of all rows of each component
    means.append(np.mean(x))

#Drag force
for j in range (len(drag_components)):
    #declaring empty array x  
    y = [] 
    
    #for i = 1; length (components) take the values of all rows in components this appends till length(components)
    y.append(df_1[drag_components[j]])
    #Averaging the values of all rows of each component
    means_1.append(np.mean(y))


#plotting 
plt.figure(1)
style.use('default')
plt.barh(components,means)
plt.title('Lift Force')
plt.xlabel('Forces(N)')
plt.grid(True, color='k')

plt.figure(2)
style.use('default')
plt.barh(drag_components,means_1)
plt.title('Drag Force')
plt.xlabel('Forces(N)')
plt.grid(True, color='k')
plt.show()

```

We would get the following plots

Lift

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Lift_1621015134.jpg)             

Drag

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/05/Drag_1621015150.jpg)   

With this information, For different Race conditions, the case setup files can be found [here](https://drive.google.com/drive/folders/1xb23ChfFNtcugdVjQ_lzY2MayNwfFco2?usp=sharing) and the conditions at which the races will be run are

Race-1

For each angle (20,45,80) the simulations will be run with a speed of 12.5m/s at the inlet and the regions will be defined with the same speed. The BC will be the law of wall for all the components of the car. The turbulence model used will be the Standard k omega 2006 model with automatic wall function. grid size will be the same as the above condition

Race-2

The simulation will be run with a speed of 20.833m/s at the inlet and the regions will be defined with the same speed. The BC will be the law of wall for all the components of the car. The turbulence model used will be the Standard k omega 2006 model with automatic wall function. grid size will be the same as the above condition                                                       

CONCLUSION:

From the above results, the following conclusions have been made.

*The underbody was offering an additional push force to the car. This result is impractical and could be the cause of poor grid size.
*The body was offering a downforce which is good. But the underbody again offers heavy lift this should practically make the car fly which is not desirable. This result could also be an error because of the change in mesh size.
*From the converge drag force plot, it is observed that at 2s and 3s there is a spike in the values of forces. This is the place where the grid changes because of grid scaling. This could indicate that the solution is grid-dependent. 
*The Pressure contour suggests that the wings are providing some downforce but the bar graph shows that this downforce is not sufficient. A change in wing design should improve the downforce. 
*The Y-plus on the body is greater than 600 which means the grid size taken is not sufficient to keep it under 100.
*Because of the less computational power and restricted cell count. It will be hard to conclude anything based on this grid size. 


 