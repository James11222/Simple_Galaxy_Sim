# Simple Galaxy Simulation
A simple illustration of how one can simulate a galaxy in python. To use this, first pip install

* `numpy`
* `matplotlib`
* `celluloid`
* `scipy`

These packages allow you to run the `galaxy.py` script. Once all necessary packages are installed, run 

`python galaxy.py`

an animated simulation of the galaxy should appear as `galaxy.gif`. On the left is the simulated positions and on the right is the simulated velocity curve.

<p align="center">
  <img src="gifs/galaxy_dark_mode.gif#gh-dark-mode-only" width="100%">
  <img src="gifs/galaxy_light_mode.gif#gh-light-mode-only" width="100%">
</p>

This simulation does not rely on derived physical models. This script is meant to demonstrate how one can simulate a galaxy in python. With careful consideration one can include more robust physical models than the ones used here. We use a characteristic constant rotational velocity curve divided by a Gaussian function centered on the Black Hole to mimic a typical velocity curve. This effectively forces the particles to behave as if they are under the influence of dark matter. We populate the galaxy with a central suppermassive black hole ~1000x more massive than the rest of the stars in the galaxy and the stars follow a decaying exponential distribution in their instantiation. Using standard Newtonian gravity we calculate the forces of gravity via a naive O(N^2) computation algorithm. This simulation is collisionless and makes use of Euler's method for computing positions and velocities.


