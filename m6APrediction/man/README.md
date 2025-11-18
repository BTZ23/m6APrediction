---
title: "README"
author: "Biting.Zhang23"
date: "2025-11-18"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# m6APrediction - m6A Site Prediction Tool

## Overview
m6APrediction is an R package for predicting m6A methylation sites in RNA sequences using machine learning models. The package provides functions for DNA sequence encoding and both single and batch prediction capabilities.

## Installation

### From GitHub
```r
# Install using devtools
devtools::install_github("BTZ23/m6APrediction")

# Or using remotes  
remotes::install_github("BTZ23/m6APrediction")
```

### From Local Source
```r
install.packages("m6APrediction_1.0.0.tar.gz", repos = NULL, type = "source")
```

## Quick Start

```r
library(m6APrediction)

# DNA sequence encoding
sequences <- c("GGACA", "ACGUU", "UACGU")
encoded_data <- dna_encoding(sequences)
print(encoded_data)

# Single prediction (requires trained model)
result <- prediction_single(
  ml_fit = your_trained_model,
  gc_content = 0.5,
  RNA_type = "mRNA", 
  RNA_region = "CDS",
  exon_length = 10,
  distance_to_junction = 8,
  evolutionary_conservation = 0.5,
  DNA_5mer = "GGACA"
)
print(result)
```

## Model Performance

The machine learning model demonstrates excellent performance in m6A site prediction:

![ROC Curve](images/roc_curve.png)
*ROC curve showing high predictive accuracy (AUC > 0.95)*

![PRC Curve](images/prc_curve.png)
*Precision-Recall curve demonstrating robust performance*

## Features

- **DNA Sequence Encoding**: Convert sequences to machine-learning features
- **Single Site Prediction**: Predict m6A modification probability
- **Batch Processing**: Efficient prediction for multiple sequences
- **Model Integration**: Compatible with Random Forest models

## License

MIT License