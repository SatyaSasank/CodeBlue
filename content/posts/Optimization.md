---
title: "Optimization of a function through genetic algorithm"
date: 2022-02-08T20:02:39+05:30

categories: [Matlab]
tags: [optimization, genetic algorithm]
toc: false
author: "Satya Sasank"
---
AIM:

The aim of this project is to write a MATLAB program to optimize the stalagmite function and find the GLOBAL MAXIMA of the function.

STALAGMITE FUNCTION: 

Stalagmite function is as shown in the figure it is just like any other function. Here we are using the stalagmite function as a function that needs to be optimized to show how we can use MATLAB’s in-built function GA to optimize the stalagmite function.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/mceclip0_1602504038.png)                                                   

                                                                   Picture credits to skill Lync

OBJECTIVES OF THE PROJECT:

Research about the genetic algorithm, optimization methods.
Plot the required graphs with titles and legends.
If possible, add any charts or workflow diagrams, or trees.

BODY OF THE CONTENT:

This project is based on Optimisation using the genetic algorithm

WHAT IS OPTIMISATION?

Optimization is the process of taking the value of a parameter to its extreme possible limits for a minimal compensation of another parameter. It is an iterative procedure. The iteration is carried until the desired value of a parameter is reached. A parameter can be something quantifiable such as length, stress, strain, force, etc. There are different methods to perform optimization such as the Deterministic method, Stochastic method. The deterministic method is used when we have all the inputs required whereas the stochastic method is used where there is a probability rather than a certain known data.

Generally, optimisation is performed in the following method

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/mceclip1_1602504075.png)                                                        

                                                               (Picture credits to Skill Lync)

The design variables are the parameters that we can vary during the optimization process.
Constraints are boundary conditions that narrow our solving by adding more functions/rules to satisfy.
The objective function is the mathematical model we build to represent the problem in a mathematical form.
Variable bounds are the range of input variables that will be varied for every iteration to get an optimized result.

Genetic Algorithm is one of the optimization algorithms. It is a search-based optimization technique based on the Evolution theory by Charles Darwin. The following diagram shows the Genetic algorithm procedure.

                                                                   
![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/mceclip2_1602504102.png)
                                                                                  Picture credits to SKILL LYNC

 The population is the number of sample data points taken to perform the optimization using the genetic algorithm. Each variation consists of huge amounts of data. In order to eliminate the unnecessary processing, we use to take samples from each variation these are called population. We could get more optimized results if we can increase the population size i.e., more samples. The given data need to be processed before usage. So, they undergo the following three steps: -

SELECTION: - In this step, the various types of data are classified according to the fitness value i.e., how well the data fits the fitness function.

WHAT IS A FITNESS FUNCTION?

A fitness function is a function that evaluates how close a given solution is to the optimum solution of the desired problem. It depends on what you define as the optimum solution. For e.g., If you were writing a genetic algorithm that simulated frog jumping, the fitness function might be the height of the jump given weight, leg size, and energy constraints. Or given a set of 5 genes, which can hold one of the binary values 0 and 1, we have to come up with the sequence having all 1s. So, we have to maximize the number of 1s as much as possible. This can be considered an optimization problem. Hence, the fitness function is considered as the number of 1s present in the genome. If there are five 1s, then it is having maximum fitness and solves our problem. If there are no 1s, then it has the minimum fitness. Using a fitness function the data is checked to see how well it fits the desired solution.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/mceclip3_1602504140.png)                                                         

                                                                                     Picture credits to skill lync

So, a Roulette wheel is made in such a way that the maximum fitting data set has more pie on the wheel. There will be a fixed point and the roulette wheel is spun to select a data set. In this way, the wheel is spun two times to get two different data sets. Usually, the more fitting data sets are selected because of their more area compared to the less fitting data set.

The selected data sets are for CROSSOVER.

CROSS OVER: - There are many ways to perform crossover. One of them is a one-point cross-over method. Refer to the image below.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/mceclip4_1602504168.png)                                                               

                                                                          Picture credits to SKILL LYNC

Here the boxes are the genes that can hold a value of 0 or 1.

MUTATION: - Mutation is nothing but a change of genes. This is carried out so that the resulting data will have diversity. If a mutation is not performed this might lead to a similar data set that didn’t result in an optimized solution and that leads to a loop. The mutation is carried out by exchanging the position of genes (0s and 1s) so that a new data set that didn’t exist with the initial data is formed.

 These 3 processes are repeated until we get the optimum value. The fitness function and the optimum criteria tell us if we have an optimum solution. Once the optimum solution is achieved the genetic algorithm is terminated. Initially, the GA progresses very fast but in the course of iterations, it tends to slow down as the improvements are very small. Hence, we terminate the GA where the solution is approximately equal to the optimum solution.

In order to use GA in MATLAB to optimize the stalagmite function, we first need to define the stalagmite function.

%This program is the stalagmite function which is mentioned in the image
%for the ease of writing the code we broke into different terms and add all
%those terms are added using the correct operator mentioned in the formula

function f = stalagmite(input_vector)

```matlab
       %since the expressions are complex to be written in a single line 
       %the expressions are broken down into small terms         
        term_11 = ((5.1)*pi*input_vector(1))+(0.5);

        f_11 = (sin(term_11)).^6;

        term_12 = ((5.1)*(pi)*input_vector(2))+(0.5);

        f_12 = (sin(term_12)).^6;

        term_21 = ((input_vector(1)-0.0667)^2)/0.64;

        term_22 = (-(4)*log(2));

        f_21 = exp((term_22)*(term_21));

        term_31 = ((input_vector(2)-0.0667)^2)/0.64;

        f_22 = exp((term_22)*(term_31));

        f = -((f_11)*(f_12)*(f_21)*(f_22)); 
       %since the GA minimizes the stalagmite function minus is added to
       %maximize the stalagmite function.
       %without the minus the graph would be improper so add do the
       %conversion before proceeding further
end
```

A genetic algorithm tends to minimize the given function. In order for us to achieve the maximum value, we need to alter the stalagmite function a little bit so that the altered function’s minimum value is the original function’s maximum value. That can be done as shown in the program above.

```matlab
clear all
close all
clc

%defining the search space
x = linspace(0,0.6,150);
y = linspace(0,0.6,150);
%defining the iterations
numcases = 60;
%Using the meshgrid in the stalgmite function results in an unknown errors
[x1, y1] = meshgrid(x, y);
%meshgrid generates a square matrix
%defining the loop
%using this loop in the stalagmite function also results in an unknown
%error so I used it here
   for i = 1:length(x1)
        for j = 1:length(y1)
        input_vector(1) = x1(i,j);
        input_vector(2) = y1(i,j);
        f(i,j) = stalagmite(input_vector);
        %calling the stalagmite function
        end
   end

% study1 statistical behaviour   
%tic and toc commands to note the time
tic
for i = 1:numcases
    %inputs denote the inputs requried for stalagmite to get a optimum of 
    %fopt, ga indicates calling the inbuilt GA function.
    %stalagmite indicates the function to be optimised
    %2 indicates the number of inputs requried for stalagmite
   [inputs, fopt(i)] = ga(@stalagmite,2);
   xopt(i) = inputs(1);
   yopt(i) = inputs(2);
end
study1_time = toc;
%the outputs at each iteration are noted and are plotted
figure(1)
%subplot command allows us to plot multiple graphs in a single window
% subplot syntax(Nth row,Mth column, Plotnumber)
subplot(2,1,1)
hold on
surfc(x1,y1,f)
shading interp% removes the grid lines formed on the graph
%plot3 command to plot 3 axis graph
plot3(xopt,yopt,fopt,'marker','o','markersize',5,'markerfacecolor','r')
title('unbounded inputs')
subplot(2,1,2)
plot(fopt)
xlabel('iterations')
ylabel('function maximum')


%study2 statistical behaviour with upper and lower bounds
tic
for i = 1:numcases
    %adding upper and lower bounds
   [inputs, fopt(i)] = ga(@stalagmite,2,[],[],[],[],[0;0],[1;1]);
   %The general syntax is 
   %[x,fval] = ga(fun,nvars,A,b,Aeq,beq,lowerbound,upperbound,nonlcon,Intcon,options)
   %where A,b,Aeq,beq= Linear Inequality constraints
   %nonlcon = non linear equality constraints
   %intcon takes integer constraints
   %where there are integer constraints GA doesnt accept equality and
   %inequality constriants
   %the terms that we dont use can be left blank []
   xopt(i) = inputs(1);
   yopt(i) = inputs(2);
end
study2_time = toc;

figure(2)
subplot(2,1,1)
hold on
surfc(x1,y1,f)
shading interp
plot3(xopt,yopt,fopt,'marker','o','markersize',5,'markerfacecolor','r')
title('bounded inputs')
subplot(2,1,2)
plot(fopt)
xlabel('iterations')
ylabel('function maximum')

%study3 statistical behaviour with upper and lower bounds along with
%population limits
%defining options
options = optimoptions('ga');
%syntax: options = optimoptions(solvername)
%defining population size
options = optimoptions(options,'populationsize',275);
%options = optimoptions(in the ga options = set the population size to 350)
tic
for i = 1:numcases
   [inputs, fopt(i)] = ga(@stalagmite,2,[],[],[],[],[0;0],[1;1],[],[],options);
   xopt(i) = inputs(1);
   yopt(i) = inputs(2);
end
study3_time = toc;

figure(3)
subplot(2,1,1)
hold on
surfc(x1,y1,f)
shading interp
plot3(xopt,yopt,fopt,'marker','o','markersize',5,'markerfacecolor','r')
title('bounded inputs with population limit')
subplot(2,1,2)
plot(fopt)
xlabel('iterations')
ylabel('function maximum')
```


The program results in different inputs required to get the maximum output each time we run the code. So, in order to solve this, the program has been divided into three parts to observe the consistency of GA with respect to constraints.

The first study is without any constraints and the time for each study is noted with the help of tic and toc commands.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/blobid0_1602504338.jpg)                                              

The second study is conducted with upper and lower bounds as we can see the maximum value region has been narrowed down by a great extent.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/blobid1_1602504353.jpg)                                              

The third study is conducted with population limit control along with upper and lower bounds and as we can see that this gives us a single solution. that matches with the surface plot.

![](https://sklc-tinymce-2021.s3.amazonaws.com/2020/10/blobid2_1602504373.jpg)                                               

The inputs in the program are manipulated to find out the correct match so that there is a minimal error in the calculation of the maximum value of the stalagmite function.

CONCLUSION:

The optimum solution of the given function is found at -1 for a given inputs of

X = 0.0668, and y = 0.0668

Reference links:

[Find minimum of function using genetic algorithm - MATLAB ga - MathWorks India](https://in.mathworks.com/help/gads/ga.html#d122e38546)


[Introduction to Genetic Algorithms](https://blog.floydhub.com/introduction-to-genetic-algorithms/)


[How to define a Fitness Function in a Genetic Algorithm? | by Vijini Mallawaarachchi | Towards Data Science](https://towardsdatascience.com/how-to-define-a-fitness-function-in-a-genetic-algorithm-be572b9ea3b4)


[How to Find Maximum and Minimum Value of a Function](https://www.onlinemath4all.com/how-to-find-maximum-and-minimum-value-of-a-function.html)
