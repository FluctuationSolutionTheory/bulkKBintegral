# README

## About
These codes are used to analyze bulk systems and calculate Kirkwood-Buff integrals (KBIs) for these systems.
```rdf.f``` calculates radial distribution functions (RDFs) and coordination numbers.
```kb.f``` calcualtes KBIs, using the coordination numbers output from the ```rdf.f``` code as input.
While below we reference a paper about electrolyle systems, these codes can also be used for non-electrolyte systems as well.
An example input file and output files are provided for the system 1M NaCl + 1M KBr (*aq*) in the ```example``` subfolder. In this example, only the (distinguishable) ion-ion RDFs and KBIs are calculated, not those involving water.
This code could also be used to calculate indistinguishable ion-ion RDFs & KBIs if the input files were setup accordingly.

## Usage
These codes require a Fortran compiler and an input trajectory file in PDB format.
The codes are setup to analyze four components. If you wish to analyze more components, change the ```ngroup``` parameter in the ```.f``` files and adjust the ```.dat``` files accordingly.

### RDF Code
To compile the code:
```sh
gfortran rdf.f -o rdf
```
Setup the ```rdf.dat``` file for the system you wish to analyze.
To run the code:
```
nohup ./rdf >& myrdf &
```

### KB Code
To compile the code:
```sh
gfortran kb.f -o kb
```
Setup the ```kb.dat``` file for the system you wish to analyze.
To run the code:
```
nohup ./kb >& mykb &
```
## Reference
E.A. Ploetz, N.D. Smyers, and P.E. Smith, "Ion-Ion Association in Bulk Mixed Electrolytes Using Global and Local Electroneutrality Constraints"
Journal of Physical Chemistry B, In Revision.

## Maintainer
[@FluctuationSolutionTheory](https://github.com/FluctuationSolutionTheory)

## License
AGPL-v3.0 © 2024 [Elizabeth A. Ploetz and Paul E. Smith]

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program (see LICENSE file).
