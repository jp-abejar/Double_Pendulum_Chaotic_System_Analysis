# Double_Pendulum_Chaotic_System_Analysis
OOP based simulation for double pendulum using parallel processing, studying poincare sections and mass variant bifurcation. Applies RK4 algorithm to estimate the differential equations of the chaotic system.

Replicating a Bifurcation Diagram for a Double Pendulum System

#Abstract
	
One very simple yet interesting system which is expected to manifest chaotic behavior is known as the double pendulum. It’s easy to see why this system is related to chaos when it is in motion since it is very difficult to predict successive positions through sole visual information. Using computational analysis, we can study this system and obtain results which are un-intuitive. One of which is a bifurcation diagram depicted in a computational physics book1 which relates the mass of the upper pendulum to the various angular velocities sampled when the lower pendulum crossed its equilibrium point. The goal here was entirely to replicate this diagram using the mechanics known to this system keeping in mind that the authors of this text are reputable and the results should be extremely similar if not exact.

#I.Background

Chaotic motion is theorized to be unpredictable or random and therefore is an interesting topic in computational analysis. In theory, given all the initial conditions for a system we should be able to accurately predict any event which solely depends on the known initial conditions. This is the case for the double pendulum which is why it is possible to accurately predict motions of this system through computational analysis. The fact of the matter is that the system is considered to be both chaotic and deterministic because of its high sensitivity to initial conditions. Computationally, we are certain of our initial conditions but for all practical purposes we can never be 100% certain of these values which ultimately leads us to label the double pendulum (among other similar mechanical systems) as chaotic. 
The most important characteristic of this system in the context in which it is studied here is the chaotic motion which it exhibits. Despite the system being considered chaotic, there are ways to analyze the data from such a system that yields peculiar patterns. Analyzing distinct properties of different trajectories can lead to bifurcations which display fractal patterns that vary depending on the variables chosen to be compared. We define fractal patterns as patterns that are recursive at infinite levels of magnification. 
The double pendulum system is composed of two pendula connected at a point which are allowed to swing freely (see Figure 1). Ideally, this system is frictionless and therefore we exclude that factors’ contribution to the motion. 

![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/dp_diagram.png?raw=true)

Fig1. A schematic of the double pendulum taken from diego.assencio.com3

Each pendulum acts as a driving force for the other so there is no need to implement driving forces in a computation. The bifurcation diagram being reproduced is seen in figure 2 where we see a pattern showing that for some mass values, the system favors a smaller set of angular velocities.



                           ![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/thtbif.png?raw=true)

Fig2. Bifurcation diagram depicting angular velocity of lower pendulum vs mass of upper pendulum.1




#II. Methods

In order to obtain the desired result, it was found that the optimal algorithm for the task would be a fourth order Runga-Kutta algorithm2. Through applying the known equations for this system, this algorithm showed to output nearly conservative total energy where the losses were less than a factor of 10-7. This algorithm is sufficiently accurate for the task at hand.  Since the text depicting the bifurcation diagram does not provide much instruction to obtain the diagram as it is displayed, manipulating the initial conditions of the system was necessary to obtain the desired outcome. The equations used for this task were obtained from a fellow classmate since the primary task was to study the output of the system rather than getting the system to function properly. The differential equations for both the energy and the motion of the double pendulum system can be seen in the fourth citation4. 
The system is initialized with zero kinetic energy, angle values where angle 1 is a small non-zero value and angle 2 is a much larger value. After testing several initial conditions, the values found to give a result closely related with the bifurcation diagram in the text were 10 degrees for angle 1 and 125 degrees for angle 2. The remaining parameters included a gravitational acceleration constant of 9.8, and both pendulums with a length of 1.0. In order to obtain various values of angular velocity for every mass, the system sampled for 70 values of angular velocity at the equilibrium point and repeated for mass values ranging from 0.1 to 14 with an interval of 0.001.  The equilibrium point here is defined where the individual pendulum angle is equal to zero (see figure 1). The angles are measured from pi to -pi, where the values are positive to the right of the equilibrium point and negative to the left of it. Since the algorithm doesn’t know our values have to be within those boundaries, it is important to set the condition to convert angle values that are beyond our boundary conditions to values that satisfy it. For this simulation, a time step of 0.001 was utilized in order to obtain better accuracy. The values for current and previous process in the algorithm were stored in two arrays for angle comparison purposes. Since the angle never exactly equaled zero, in order to obtain an angular velocity at that point, a condition was set to identify where the previous angle and the current angle were on opposite sides of the equilibrium point. When this was the case, the angular velocity at the current point was saved to an array which would later be graphically compared to an array holding all the values set for the mass range. The resulting image depicted a bifurcation diagram. 











#III. Results




![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/bifmult.png?raw=true)

Fig 3. Bifurcation diagram attempts (a-f) where the initial angles were varied to analyze and compare output to the desired diagram. Initial angle 1 was modified in every instance keeping all other initial conditions constant. 


![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/bif.png?raw=true)

Fig 4. Bifurcation diagram (partially enhanced) for double pendulum displaying the relationship between angular velocity of the lower pendulum as it crossed the equilibrium point 70 times and the mass of the upper pendulum. Mass values ranged from 0.1 to 14 with an interval of 0.01 to display with greater detail.

![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/bif2.png?raw=true)

Fig 5. Bifurcation diagram for double pendulum focused at the beginning mass values to show detail.


![alt text](https://github.com/jp-abejar/Double_Pendulum_Chaotic_System_Analysis/blob/main/Img/bif3.png?raw=true)

Fig 6. Bifurcation diagram for double pendulum focused at the center mass values to show detail.

The graphs displayed in figure 3 were plotted to determine what initial conditions were necessary to obtain the desired relationship between the angular velocity of the lower pendulum and the mass of the upper pendulum. The mass intervals were set at 0.5 in order to get a general relationship to decide what initial values would be given to the system to process for a smaller interval of mass values. We can see that as we increase angle 1, our diagram becomes more and more similar to the one we want to replicate, but as it increases past the angle of 12.5 degrees, the diagram no longer tends to resemble the diagram in figure 2. Figures 4 through 6 depict the results of a bifurcation diagram shown at full scale, and focused in at different sections to show pattern details. The only initial conditions manipulated were the starting angle positions which were set at 10 degrees for angle 1 and 125 degrees for angle 2. 



#IV.Summary

The results obtained in this experiment did not exactly match the desired diagram but did show very similar qualities to it. The text diagram shows a graph that converges toward the approximate value of 6 rad/sec as the mass of the upper pendulum is increased, which matches the obtained diagram. We also see that for both graphs, the system tends to periodically favor a small set of velocities as the mass increases. This pattern of data densities is what leads to the conclusion that the bifurcation diagram displayed in the text is very well plausible given the proper initial conditions. Since the system is extremely sensitive to initial conditions, we can expect other similar bifurcation diagrams for other initial conditions, but we can clearly see from the diagrams in figure 3 that the convergence to a small range of velocities with the increase of upper pendulum mass is not possible for all combinations of initial conditions since there appears to be a limit on how the initial angle values can be related to each other to obtain this type of graph. The specific conditions for the initial values that output a diagram like that of figure 2 are undetermined for this experiment. It is only seen that a plausible relationship between the initial starting angles that gives a fairly similar result is that the angle of the first pendulum should be close to its equilibrium point, while the second angle should start at a large angle.





1 Landau, Rubin H., et al. Computational Physics: Problem Solving with Python. (pp.377) ,Wiley-VCH, 2015.

2 Cheever, Erik, and Swarthmore College. “Fourth Order Runge-Kutta.” The Laplace Transform, lpsa.swarthmore.edu/NumInt/NumIntFourth.html.

3 Assencio, Diego. “The Double Pendulum: Lagrangian Formulation.” Funnel: How Long Does the Fluid Take to Go through? - Diego Assencio, diego.assencio.com/?index=1500c66ae7ab27bb0106467c68feebc6.


4 “Double Pendulum.” Maritime Theater, web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html.

