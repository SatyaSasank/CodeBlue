---
title: "Flow over ahmed body"
date: 2022-02-08T16:08:18+05:30

categories: [Ansys]
tags: [CFD, Aerodynamics]
toc: false
author: "Satya Sasank"
---
AIM:

This project aims to simulate the flow over Ahmed body in Ansys Fluent.

OBJECTIVES OF THE PROJECT:

*Generate Velocity and pressure contours.  

*Calculate the drag coefficient plot for a refined case.  

*Generate the vector plot clearly showing the wake region.   

*Perform the grid independency test and evaluate the values of drag and lift with each case.    

THEORY

What is an Ahmed body?

The Ahmed body is a geometric shape first proposed by S. R. Ahmed and Ramm in 1984. The shape provides a model to study geometric effects on the wakes of ground vehicles (like cars). It is used to validate the CFD solvers. The results of airflow over this body were captured using the wind tunnel experiments and the CFD solver results are compared with these experimental results to find out the accuracy of the solver.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/AhmedBody_1623386070.png)

Picture credits : [Flow and Turbulence Structures in the Wake of a Simplified Car Model](https://www.researchgate.net/publication/266883948_Flow_and_Turbulence_Structures_in_the_Wake_of_a_Simplified_Car_Model_Ahmed_Model)


The above picture represents the ahmed body dimensions. Its virtual version can be found [here](https://drive.google.com/open?id=1gfeLjCjeeVqio38CGdhDEsfRUotrU5l_)


PROCEDURE:

Importing the CAD model

The geometry can be created using any standard CAD package. In this case the geometry is already created.

We import the STEP file into Ansys Spaceclaim. Once the geometry appears, we need to check the geometry units. We can do this by selecting any one side of the geometry. Upon selecting the dimensions are displayed at the bottom. If the dimensions are in mm.

We can change the dimensions using File> Spaceclaim Options > Units. And convert them to meters.

To observe the flow over an Ahmed body we create a virtual wind tunnel. We do this by creating another geometry on this. This can be created by using enclosure.

But before we create the wind tunnel, we need to know the dimensions of the wind tunnel.

The distance after the Ahmed body should be 5 - 10 times the characteristic length but 10 times the characteristic length requires more computational time and power hence we reduce it to 5 times to accurately capture the wake produced.

In the z-direction, the total length should be at least 5 times but due to computational restraint we reduce it to 0.5 m

we keep the same length for y as well, Y direction length (from above the Ahmed body) is 1 m

Now we use the enclosure option to create an enclosure on the ahmed body. We then uncheck the original imported geometry, suppress it for physics. This is to be done to prevent overlapping of bodies.

Now we create another box around the ahmed body, This box will be used for mesh refinement. We use the interference tool to remove the mesh interference. Once this is done we exit the Spaceclaim. This model will be updated to Ansys workbench.

Once created the geometry should like this

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(257)_1623415152.png)

Meshing

Now we open the Ansys mesher. Here we first start by naming the boundaries as.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Boundary_name1_1623415495.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Boundary_name2_1623415534.png)         

The second box created is for refining the mesh hence we select the mutizone mesh. If we inspect the geometry we can see that the legs of the ahmed body lost their cylindrical shape. This is due to coarse size of the mesh. To regain the shape we use face sizing at the legs. In order to give the refinement we use the edge sizing option to refine the mesh inside the refinement zone. To capture the effects at the car wall we give inflation layers at the car wall. Upon meshing the geometry will look like this.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(260)_1623478493.png)

Case setup


For this project we use a steady state solver. And since the velocity is less to be worried about mach number we use a pressure based solver. This is an external flow hence we use a K Omega SST model to solve for turbulence. The fluid is taken as air and the properties are default. The boundaries are set as follows.

![](images/ahmedbodyTable.png)

Before we proceed to initialisation we need to take care of the reference values. This is important and incorrect reference values will lead to improper drag and lift values.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/ReferenceValues_1623480486.png)


Once these values are set we can proceed to initialisation (Hybrid/Standard). Once initialisation is done we get initial velocity contours. We create a plane and observe the velocity contour on that plane. Also we save this animation for future references.

We now run the simulation for 350 iterations and the end results are as follows

Velocity contour
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Velocity_contour_1623484042.png)

Residual plots
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Residuals_1623484066.png)

Cd plot
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Cdplot_1623484109.png)

Cl plot
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Clplot_1623568874.png)

POST PROCESSING:

We can post process the above results to get some important results such as line probes, graphs and vector plots. These are important to identify the point of seperation, recirculation zone etc.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/VectorPlot_1623498785.jpg)

The above vectors are plotted on a plane XY at Z = 163.5 mm.

Recirculation Zone

When an object moves in a fluid, the first layer of the fluid sticks on the surface of the body. When there is a sharp change in the geometry, 
there comes a point on the body where the fluid detaches from the surface, thus breaking the boundary layer, this point at which the fluid 
detaches from the surface is called the point of separation. Due to this flow seperation happens and a low pressure zone is created at the 
end of the ahmed body.This low pressure is again occupied by air from the high pressure. This process is called as recirculation of air and this
zone is called as ReCircualtion Zone 
If we are to take a look at the pressure contour it would look like this. 

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(269)_1623495759.png)

Negative Pressure ?

In the above contour, the range of pressure is adjusted so that the negative pressure range can be clearly seen. The explanation as to why there is a negative pressure in the domain goes as follows:

When the air hits the ahmed body the air gets stagnated at the front causing a stagnation of air due to which the pressure rises, This can be observed at the front portion of ahmedbody. The area in red indicates the increased pressure.
This air eventually passes over the geometry and when ever there is a sharp change in geometry this air tends to expand increasing velocity and losing pressure. This pressure drop causes the negative pressure zones.

Grid Dependency Test

The above case we looked at is the most refined case. We change the mesh size and see if it affects the solution. In the earlier case the size of the grid in the refinement zone is 0.02m. We change the size of the refinement zone to 0.03 and 0.04 m

GRID SIZE 0.03 m
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(277)_1623567293.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(278)_1623567304.png)

GRID SIZE 0.04 m

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(273)_1623566405.png)
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/06/Screenshot%20(279)_1623568424.png)


CONCLUSION:

*From the Velocity contour, it is observed that after passing through the geometry there is a region where recirculation happens this region is called as wake region. Due to this low pressure, the air at the front tries to come back which causes the Pressure drag.  
*On reducing the grid size the Cd and Cl values increases by some extent implying the solution is grid dependent. But the minor increase indicates that by the next lower size the solution may stabilise and become grid dependent.  

REFERENCES:

[The Ahmed Body](https://curiosityfluids.com/2016/09/07/the-ahmed-body/)

[The viscosity of Air, Dynamic and Kinematic](https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm#:~:text=The%20viscosity%20of%20air%20depends,10%2D5%20Pa%C2%B7s%20.&text=At%2025%20%C2%B0C%2C%20the,the%20kinematic%20viscosity%2015.7%20cSt.)

[Experiments and numerical simulations on the aerodynamics of the Ahmed body](https://www.researchgate.net/publication/330383775_Experiments_and_numerical_simulations_on_the_aerodynamics_of_the_Ahmed_body)




