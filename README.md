# droplet-microfluidics-lbm
Modeling of droplet dynamics using multi-phase lattice Boltzmann method

This project presents the coalescence and splitting dynamics of droplets for two-phase fluids. The
Shan-Chen model (10.1103/PhysRevE.47.1815) is used to capture the interaction between the different phases, including the transfer of
momentum, energy, and mass.

The model is tested for an exemplary microfluidic channel, with rigid and stationary walls (bounce-back boundary condition) and periodicity in x-direction.

The animation below shows the transport of the droplets inside the microfluidic channel to perform the first coalescence and then splitting of droplets;
and the color bar denotes velocity of the base fluid, and therefore maximum value is present at the interface. For the simulation details, two droplets are placed at the beginning of the channels. For t < 1300, a flow $\Delta p > 0$ is introduced to the fluid at the entrances of two channels to allow the droplets to meet at the junction. As the coalescence is observed, the flow is resetted $\Delta p = 0$ to let the larger droplet undergo relaxation, for 1300 < t < 1800. After t > 1800, a reversed-flow $\Delta p < 0$ is introduced to demonstrate splitting the larger droplet into two parts.

<p align="center">
    <img src="[channel.gif](https://github.com/lynspica/droplet-microfluidics-lbm/blob/main/figs/channel.gif)" alt="animated" />

![](https://github.com/lynspica/droplet-microfluidics-lbm/blob/main/figs/channel.gif)
</p>
