# Simultaneously fitting equilibrium and dynamic data for diffusion-controlled sorption

The Matlab codes included herein are used for estimating the effective diffusivity of diffusiving species for indivual step-changes in relative pressure, and then simulataneously fitting solution diffusion + immobilizaiton model to the effective diffusivity and equilibrium isotherm data simultaneously. The code is designed to accept user-selected inputs of IGAsorp dynamic sorption datafiles but can be manually supplied arrays of time-mass-activity-temperature data as well.

## Getting Started

Download DeffCalculator.m and FitIsoAndDeff.m to a local Matlab folder directory.

### Prerequisites

Matlab R2020a or newer. These can be obtained from mathworks.com. Alternatively, users without access to Matlab may consider using GNU Octave, a free open-source clone of Matlab, available at https://octave.org. For use in GNU Octave, the code will have to be modified accordingly by the user. There is no guarantee of retaining all code functionality when porting the code for GNU Octave.

## Getting Started Example

Download Raw_Data_Compiled_Excel.xlsx from the supporting information of "A general model form for capturing concentration-dependent effective diffusivities in a range of composite polymeric materials." Open the Excel file, and on the "Table Of Contents" tab click the hyperlink for Polysulfone 20 °C Dynamic Curves. Select cells A1:H7012 and copy to a .txt file in your preferred folder. In Matlab, ensure DeffCalculator.m is in the Matlab working directory and run this function by typing "DeffCalculator([0:20:80,90])" into the command line. The code will prompt the user to supply Sample Name (Polysulfone), the sample geomtry (0 for slab), and the sample radius/thickness in cm (0.176). A dialog window will open requesting the use to "Select the Input Dynamic File(s)." Navigate to the saved text file and select "Open." The code will then extract the data from the text file and parse each humidity step based on the user inputted array. In this example, an input array [0:20:80,90] informs the code to parse the dynamic data at 20%, 40%, 60%, 80%, and 90% during sorption and desorption.

The code will create figures of each humidity step and show the progress of fitting a effective diffusivity to the experimental data. Once completed, the effective diffusivity data will be saved in a text file names "[filename]_Diffusivity.txt" and should look like the example below.

Polysulfone

Thickness/Radius: 0.176 cm

Geometry:Slab

Temperature (°C):20

RH	D / cm^2/s	Temperature(°C)	Weights

0.200000	3.859221e-08	20.000000	1.000000

0.400000	3.610725e-08	20.000000	1.000000

0.600000	3.862852e-08	20.000000	1.000000

0.800000	3.811194e-08	20.000000	1.000000

0.900000	4.137525e-08	20.000000	1.000000

0.800000	4.252364e-08	20.000000	1.000000

0.600000	3.036380e-08	20.000000	1.000000

0.400000	3.727293e-08	20.000000	1.000000

0.200000	3.745932e-08	20.000000	1.000000

0.000000	3.399344e-08	20.000000	1.000000

Navigate to the "Polysulfone isotherms" tab in the Excel file (manually or through the hyperlink in the "Table of Contents") and copy cells A1:D30 to .txt file in a desired folder. Run FitIsoAndDeff.m from Matlab by typing "FitIsoAndDeff" in the command window. The code will prompt the user for the sample name (Polysulfone) and then request the user to supply which sorption modes should be included in the mobile or immobile concentrations. For this example, when the first prompt appears:

"MOBILE sorption modes: Langmuir (1/0) Henrys (1/0) Pooling (1/0) entered as [L,H,P]:"

we enter an array [0,1,0] which tells the code we only want mobile Henry's mode for polysulfone. Following the second prompt:

"IMMOBILE sorption modes: Langmuir (1/0) Henrys (1/0) Pooling (1/0) entered as [L,H,P]:"

we enter an array [0,0,0] to indicate that we do not want any immobile species for this material. In general, the user can test multiple models and adjust which sorption modes are present in mobile or immobile concentrations as necessary for the material.

A dialog box will now open first requesting users to select the isotherm files (hold control to select multiple files at once, if desired) and select "open" to import these files into Matlab. Next another dialog box opens to request users to select the diffusivity files (which will end in "_diffusivity.txt"). The code will now employ the Latin Hypercube method to sample the variable space with 100 initial guesses, and minimize each guess to fit the model to the isotherm and effective diffusivity data simulataneously. Each time a better solution is found, the code prints to the command window the iteration number and the current best weighted L2 error. The final best fit solution is saved in a text file in the same directory as the isotherm files with the name "Iso+Dyn Parameters.txt". For this example, which used only the 20 °C data, the paramters file will show Hm(T=20 °C)=8.02 mg/cm3 and D(T=20 °C)=3.74e-08 cm2/s. With more data supplied, the code will also provide H and D as functions of temperature and the 95% t-confidence intervals for each of the parameters.

## What the codes do

DeffIso uses the matlab built-in pdepe function to solve the 1-D diffusion equation da/dt = x^-m d/dx Deff\*x^m da/dx where x is the spatial variable and m = 0 for slab, 1 for tall cylinders, and 2 for spheres. Note that thin cylinders are treated as slab. The interior boundary condition at the symmetry line of the material is da/dx = 0 and at the boundary is a=a(t), where the real relative humidity or concentration data is used as the boundary condition and is non-constant. This non-constant boundary condition can slow down the matlab code becuase it will evaluate the boundary condition more often if higher time-resoltuion data is supplied. The time resolution supplied to the code can be altered by defining the samplerate variable, for example, DeffCalculator([0:20:80,90],10), will only include ever 10th datapoint when fitting Deff. This can greatly speed up the estimatation of Deff when time resolutions are excessively high.

The code estimates Deff by first guessing Deff based on the slope of the uptake versus sqrt(time) curve, and then finds lower and upper bounds on the best fit. The bisects the lower and upper bound iteratively until the certainty in the best-fit Deff is below the desired threshold (maxerror). The user can choose to supply a Matlab array as the dynamic data, by 

DeffCalculator(RH_steps,[],[],data)

where empty arrays skip assigning values to samplingrate and maxerror. The data array must have four columns in the following order:time (in minutes), mass, activity, temperature (°C).

FitIsoAndDeff simultaneously fits the c(a,T) and Deff(a,T) using the Latin Hypercube method to generate 100 initial guesses and minimizing each initial guess to find the minimum weighted L2 error. The error for each data point is normalized by a rough estimate of their relative variance, which is taken as the square of the median for each dataset (concentration or Deff) and concentrations are weighted by an additional factor of 100 to recognize the higher confidence we have in their precise values over the effective diffusivities.

The model fits to the data are in the form:

total concentration = mobile(Langmuir + Henry's + Pooling) + immobile (Lanmguir + Henry's + Pooling)
effective diffusivity = mobile concentration/a/(d(total concentration)/da)\*D

where D is the intrinsic diffusivity and a is the activity. Activity can be relative pressure, pressure, or concentration. Further:

Langmuir = b\*a/(1+b\*a)\*L

Henry's = H\*a

Pooling = p\*a^n

with separate variables used for mobile and immobile species. These codes can fit multiple temperatures simultaneously if the user supplies them (by holding control when selecting isotherm and diffusivity files), in which case the following variables take Arrhenius or van't Hoff forms: b, H, p, D. These temperature dependencies are described by equations of the form:

f(T) = f(T = 20 °C)\*exp(-dlnfdinvT\*(1/T - 1/293.15 K ))

where f is a variable (b, H, p, or D) and dlnfdinvT is the derivative d ln f/d(1/T), which is related to the activation energy (or enthalpy of sorption) by Ea (or dH) = -R d ln f/d(1/T). When fitting multiple tempratures, f(T = 20 °C) and dlnfdinvT are two independent parameters.

Confidence intervals for parameters are determined by calculating the jacobian (J), approximating the Hessian (H) as H=inv(transpose(J)\*J). The standard error squared is the error / dof, where dof is the degrees of freedom = # of datapoints - # of parameters. The t-value corresponding to the 95% confidence interval is an inverse of the Student-T function for teh specified degrees of freedom at a value of 0.025. Finally, the error for each parameter is t_value\*sqrt(standard error squared \* diagonal(H)), where each diagonal of the Hessian corresponds to variance of a parameter.

## License
MIT License

Copyright (c) 2024, Lawrence Livermore National Security, LLC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Release

LLNL-CODE-863665
