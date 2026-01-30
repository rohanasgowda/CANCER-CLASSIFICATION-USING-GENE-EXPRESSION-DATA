# Cancer Classification Using Gene Expression Data

[![MATLAB](https://img.shields.io/badge/Language-MATLAB-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Status](https://img.shields.io/badge/Status-Completed-success)](#)

## 📌 Overview
This project implements a **binary cancer classification system** using gene expression data from NCBI (GSE9476). The goal is to classify samples as **cancerous (AML)** or **control**, using machine learning with **feature selection** and **hybrid optimization** to improve accuracy.

![Feature Selection Funnel](./Feature%20Selection%20Funnel.png)

## 📂 Repository Structure

| File | Description |
| :--- | :--- |
| `cancerclassifier.m` | **Core Script:** MATLAB implementation of the classification pipeline. |
| `GSE9476_series_matrix.txt` | **Dataset:** Cleaned gene expression profiles from NCBI GEO. |
| `LICENSE` | **Apache License 2.0:** Permissive license with explicit patent grants. |
| `Confusion Matrices.png` | Performance breakdown for each model variant. |
| `Feature Selection Funnel.png` | Visualization of the gene reduction process (12k → 200 → 87). |
| `Model Accuracy Comparison.png` | Bar chart comparing SVM, mRMR, and RF+SMO Hybrid results. |
| `PCA of SMO-Selected Genes.png` | Scatterplot showing sample clustering based on optimal features. |

---

## 🧩 Tools & Libraries
* **MATLAB**
* **Machine Learning Models:**
    * Support Vector Machine (SVM)
    * Random Forest (RF)
* **Feature Selection Techniques:**
    * mRMR (Minimum Redundancy Maximum Relevance)
    * Spider Monkey Optimization (SMO)

![Model Accuracy Comparison](./Model%20Accuracy%20Comparison.png)

## ⚙️ Workflow
1.  **Data Loading:** GEO dataset (GSE9476) is loaded and cleaned.
2.  **Hold-Out Split:** 70% training / 30% test.
3.  **Model 1 – SVM:** Using all genes (~12,625).
4.  **Model 2 – SVM + mRMR:** Select top 200 genes using mRMR.
5.  **Model 3 – RF + SMO Hybrid:** Random Forest ranks top 200 genes; SMO selects optimal subset (~87 genes).
6.  **Evaluation:** Accuracy, confusion matrices, and visualizations.

![PCA of SMO-Selected Genes](./PCA%20of%20SMO-Selected%20Genes.png)

## ✅ Results
* **SVM (All Genes):** ~78.95% accuracy
* **SVM + mRMR (Top 200):** ~84.21% accuracy
* **RF + SMO Hybrid:** **~94.74% accuracy** with 87 genes selected

![Confusion Matrices](./Confusion%20Matrices.png)

---
**Developed by [rohanasgowda](https://github.com/rohanasgowda)**
