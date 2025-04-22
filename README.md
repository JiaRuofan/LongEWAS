# Overview

This repository contains the source c code for the paper Detection of Cell-type-specific Differentially Methylated Regions
in Epigenome-Wide Association Studies. The code will later be packed into a R package

# Introduction

DNA methylation at cytosine-phosphate-guanine (CpG) sites is one of the most important types of epigenetic mechanisms and has been associated with many phenotypes. As DNA methylation levels change substantially throughout life due to development and aging, longitudinal epigenome-wide association studies (EWAS) are emerging. However, despite the active research on detecting cell-type-specific associations for cross-sectional EWAS data, rigorous statistical methods for analyzing longitudinal EWAS data are lacking. As EWAS data measures the DNA methylation levels at the bulk level instead of the single-cell level, the obtained methylome for each sample is the signals aggregated from different cell types. Therefore, in order to model the cell-type-specific methylation trajectories, random effects have to be added to the unobserved latent cell-type-specific methylome that needs to be deconvoluted from the aggregated-level observed data instead of being added to the outcome layer as common longitudinal studies do, which makes the statistical inference challenging. Here, we propose a hierarchical model, LongEWAS, to perform cell-type-specific association detection for longitudinal EWAS data. We prove the model identifiability for LongEWAS, develop a generalized expectation-maximization algorithm for parameter estimation and use the Louis estimator to construct test statistics. 

# Maintainer

Ruofan Jia 1155165249@link.cuhk.edu.hk
