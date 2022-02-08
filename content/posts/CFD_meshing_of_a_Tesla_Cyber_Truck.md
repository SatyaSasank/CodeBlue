---
title: "CFD meshing of a Tesla Cyber truck"
date: 2022-02-08T19:32:12+05:30


categories: [ANSA]
tags: [Pre-Processing, Meshing]
toc: false
author: "Satya Sasank"
---
AIM:

The aim of this project is to pre-process the tesla cyber truck and mesh it. The CAD model of which can be found [here](https://drive.google.com/file/d/1EqotfJB3Ur_WdTbkW1J3DUns-ngAWqg4/view?usp=sharing).

OBJECTIVES:

For the given model, check and solve all geometrical errors and Assign appropriate PIDs. Perform meshing with the suitable Target length and element Quality criteria.

PROCEDURE:

We start with importing the file.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(615)_1632849475.png)

The geometry is filled with single cons to eliminate these we use the Topo option.

Before we get into TOPO cleaning procedure, we have to slice the geometry into two parts in their symmetrical plane, since they are symmetrical in nature also.

On completing all the procedures in one half of the model, we can make a copy of this half into the other easily at last.

Use the delete option to delete one-half of the model.

As we see the model given has many single cons, we want to do some Geometry checks to explore what are the errors we are encountering here,

After completed with the TOPO clean up the geometry is shown below,

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(616)_1632850080.png) 

After this step, we are flagging the appropriate boundaries to the respective regions as shown below,

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(617)_1632850356.png)

After giving the boundaries we are meshing the geometry with the given mesh lengths.

The meshed Geometry is shown in the below picture,

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(618)_1632851130.png)

Once the mesh is generated on all the regions of the truck, we can now generate the complete truck model.

This is done by the use of the symmetry option. Transform > copy > entities > symmetry > mirror 3 point plane, then select any 3 points on the symmetric plane.

Once this is done, all the entities will get copied to the other side of the symmetric plane.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(619)_1632851264.png)

Once this is done, the truck geometry must be enclosed in a box (virtual wind tunnel) to perform external aerodynamic simulations. If the length of the car is 'L' meters then, the virtual wind tunnel dimensions are as follows:

a. Distance from the wind tunnel inlet to the front of the truck= 4L

b. Distance from the aft of the truck to the wind tunnel outlet= 6L

c. Distance from the sides of the truck to its respective wind tunnel sides= L

d. Distance from the floor to the top surface of the wind tunnel= 3L

The virtual wind tunnel is created by first, identifying the origin of the global coordinate system, then, the location of wind tunnel floor corners is calculated and points created, these points are joined with curves.

A face can be created from these curves. To get the top surface of the wind tunnel, a copy of this face is offset at a distance of 3L.

The side faces and the inlet and the outlet of the wind tunnel can be created by using the reference of the wind tunnel top surface and floor cons.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(620)_1632851367.png)

Now a surface mesh must be generated on the wind tunnel surface.

Since we want finer volume mesh near the surface of the car, which is resting on the floor of the wind tunnel, we must provide a smaller mesh size for the wind tunnel floor and larger mesh size on the wind tunnel top surface.

The mesh size on the sidewalls and the inlet and the outlet of the wind tunnel must increase gradually from the floor to the top surface of the wind tunnel.

To achieve this first the wind tunnel floor and top surfaces are meshed separately with triad mesh elements and the mesh target length of 20mm and 100mm respectively.

Now, to mesh the other faces in a gradually increasing manner, go to mesh> perimeters> spacing, select the 4 edges, take biasing as geometric (1.1), and dmin= 20, dlimit= 100.

Next, mesh the geometry using spot mesh, and check for quality criteria errors if any. If there are any errors, it can be rectified mesh> shell mesh> reconstruct.

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(621)_1632851535.png)

![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/09/Screenshot%20(622)_1632851547.png)

Once the surface mesh is generated on all the components, the volumetric mesh can be created.

In volume mesh> define> auto-detect> whole data base.

Here, delete all the unwanted components, ie, the volumes enclosed by the solid surfaces.

Now right click on the fluid volume detected and click on re-mesh using tetra CFD mesh, and the required mesh will be generated.

CONCLUSION:

The tesla cyber truck model has been meshed with the given criteria. A total of 35792848 2D elements have been generated. This is a huge number of elements and is not ideal in most cases hence the mesh size will be adjusted to get the optimal trade-off between the mesh size and the number of elements. The generation of 3Delements for this model requires huge power resources and time. Hence the volumetric meshing of this model cant be done in this project.  

The final result of this project can be found [here](https://drive.google.com/file/d/1yc0hhtO2GTgD0c8AIcJjGT0jWUkU6kNM/view?usp=sharing).