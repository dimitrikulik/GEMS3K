# GEMS3K

Numerical kernel solver of the GEM-Selektor v.3 package for geochemical modelling. 
Just extended with an alternative Activity C++ API as the basis for a new Python API (based on Pybind11) on top of it, see xGEMS repository.

### Briefly about GEMS3K

* Implements the improved GEM IPM-3 algorithm with excellent mass balance precision and fast convergence to Gibbs energy minimum even in very complex non-ideal chemical systems with two-sided metastability constraints (learn more on GEMS3K web page).
* Written in C/C++. Using compiler directives, the GEMS3K code can be compiled as a standalone program e.g. 'gemcalc'; as a static or dynamic library for coupling with the mass transport simulator or another code; or as part of the GEM-Selektor v3 code together with GUI and databases.
* Version: currently 3.4.1.

### License
* GEMS3K is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* GEMS3K is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
* You should have received a copy of the GNU Lesser General Public License along with GEMS3K code. If not, see http://www.gnu.org/licenses/. 

### How to clone (download) GEMS3K source code

* In your home directory, make a folder named e.g. ~/gitGEMS3 with a subfolder gitGEMS3/standalone.
* Change into ~/gitGEMS3/standalone and clone this repository from https://bitbucket.org/gems4/gems3k.git using a preinstalled free git client SourceTree or SmartGit (the best way on Windows). 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line (do not forget a period) to download the actual "trunk" branch:
~~~
git clone https://bitbucket.org/gems4/gems3k.git . 
~~~
* To switch to another branch (e.g. devEJDB), use a git client or open a terminal, cd to ~/gitGEMS3/standalone, and type in the command line
~~~
git checkout -b branches/devEJDB --track origin/branches/devEJDB
git pull origin branches/devEJDB
~~~
To switch back to trunk, type
~~~
git checkout trunk
~~~

### How to build GEMS3K library and examples ###

For details on how to build the GEMS3K chemical solver library and examples, please consult this web page:

http://gems.web.psi.ch/GEMS3K/techinfo.html

TBD
* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines