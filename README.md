NB Github markdown does not seem to natively support equations. That's why they look weird. \[TODO: Fix this\]


# PWP model

File: PWP.py

Usage: PWP_tests.ipynb

This is code for the PWP ocean mixed layer model (Price, Weller, and Pinkel, 1986). The main function is model_timestep, which takes the following inputs:
* T: Temperature profile (deg C)
* S: Salinity profile (ppt)
* U: Zonal velocity profile (m s^-1)
* V: Meridional velocity profile (m s^-1)
* z: Depth profile (m; must be evenly spaced) TODO: allow this to be not evenly spaced
* I: Short-wave radiation (W m^-2, positive incoming), which is taken up by the ocean as $$I1 e^{-z/lambda1} + I2 e^{-z/lambda2}$$
* L: Long-wave radiation (W m^-2, positive incoming), which is taken up/emitted by the uppermost layer
* E: Evaporation (kg m^-2 s^-1), positive-definite), which is treated as added salinity to the uppermost layer
* P: Precipitation (kg m^-2 s^-1),positive-definite), which is treated as removing salinity from the uppermost layer
* tau_x: Zonal wind stress (N m^-2), spread evenly through the mixed layer
* tau_y: Meridional wind stress (N m^-2), spread evenly through the mixed layer
* dt: Timestep (s)
As well as the optional inputs: \[TODO: re-organize these in a better way\]
* use_static_stability (boolean, default=True): if true, PWP model ensures a statically stable water column (see below)
* use_mixed_layer_stability (boolean, default=True): if true, PWP model ensures a stable mixed layer (see below)
* use_shear_stability (boolean, default=True): if true, PWP model ensures a water column with stable vertical shear (see below)
* use_Ekman_flux (boolean, default=False): if true, PWP model adds in a user-defined Ekman heat flux
* use_MLI (boolean, default=False): if true, PWP model adds in a user-defined mixed layer instability restratification heat flux
* Ekman_Q_flux (W m^-2, default=None): Ekman forcing parameterized submesoscale heat flux (see below), only used if use_Ekman_flux is True
* MLI_flux: (W m^-4, default=None): When multiplied by the square of the mixed layer depth, this is the parameterized mixed layer instability heat flux (see below), only used if use_MLI is True
* verbose (boolean, default=False): if true, model will print out a lot of annoying numbers, but could be useful to debugging
* return_MLD (boolean, default=False): if true, model additionally returns the mixed layer depth
Finally, the model needs a bunch of other parameters to run. These all have defaults, but can be defined as well:
* I1 (unitless, default=0.62) fraction of shortwave radiation absorbed with a decay scale of lambda1
* I2 (unitless, default=None) fraction of shortwave radiation absorbed with a decay scale of lambda2; if not set, will revert to 1-I1
* lambda1 (m, default=0.6) decay scale of I1 fraction of shortwave radiation
* lambda2 (m, default=20) decay scale of I2 fraction of shortwave radiation
* T0 (deg C, default=17) reference temperature
* S0 (ppt, default = 36) reference salinity
* rho0 (kg m^-3, default=None) reference density, if not set will be calculated by T0 and S0
* alpha (kg m^-3 degC^-1) thermal expansion coefficent (should be negative), if not set will be calculated by T0 and S0
* beta (kg m^-3 ppt^-1) saline contraction coefficient (should be positive), if not set will be calculated by T0 and S0
* f (s^-1) Coriolis frequency, if not set will be calculated as the frequency at 40 deg N \TODO: make sure okay if this is in the S.H.\]

_Stability types_ \[TODO: Make this more clear\]

*Static stability*: Ensures that $$\rho_z \geq 0$$. If this is not fulfilled, the water column is mixed to deeper depths until this is true.

*Mixed layer stability*: Ensures that $$R_b \geq 0.65$$, where $$R_b = \frac{g\delta\rho h}{\rho_0(\delta U^2 + \delta V^2)}$$ is the bulk Richardson number. If this is not fulfilled, the water column is mixed to deeper depths until this is true.

*Shear stability*: Ensures that $$R_g \geq 0.25$$, where $$R_g = \frac{g\rho_z}{\rho_0 (U_z^2+V_z^2)}$$ is the gradient Richardson number. If this is not fulfilled, water masses are linearly mixed together over increasing depths until this is true.

_Submesoscale heat fluxes_ \[TODO: Make this more clear\]

*Ekman heat flux*: Thie is calculated as $$\frac{\tau\times\nabla b}{f}\frac{c \rho_0}{|\alpha| g}$$

*MLI heat flux*: This is calculated as $$C\frac{\nabla_b^2 h^2}{|f|}\frac{c \rho_0^2}{|\alpha| g}}$$
