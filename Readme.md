## Table of contents
* [General info](#general-info)
* [Models](#models)

## General info
This project implements Baum-Welch algorithm to predict trees from output file of Mixtures of trees(MAT) model.
The input file for B-W algorithm is .sitelh output file from iqtree.
Project is created using R, depends on R package "aphid" and "testit".
The output is vector indicating the tree for a given site.

## Models
Model 1 - The tree with highest probability at a given site is considered, and converted into a sequence for training of B-W algorithm.
Model 2 - The tree with highest probability at a given site and constant sites are considered, and converted into a sequence for training of B-W algorithm.
Model 3 - The tree with highest probability at a given site constant sites and non-informative sites are considered, and converted into a sequence for training of B-W algorithm.
Model 4 - The tree with highest probability at a given site constant sites, non-informative sites and same parsimony sites are considered, and converted into a sequence for training of B-W algorithm.
Mix Model - Starts from model 1 and moves to next model if B-W algorithm doesn't converge in 100 iterations.
Gen-model - Model 4 generalized for n number of trees.
