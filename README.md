# Cancer Classification Using Gene Expression Data

## 📌 Overview
This project implements a **binary cancer classification system** using gene expression data from NCBI (GSE9476). The goal is to classify samples as **cancerous (AML)** or **control**, using machine learning with **feature selection** and **hybrid optimization** to improve accuracy.

---

## 🧩 Tools & Libraries
- MATLAB  
- Machine Learning Models:
  - Support Vector Machine (SVM)  
  - Random Forest (RF)  
- Feature Selection Techniques:
  - mRMR (Minimum Redundancy Maximum Relevance)  
  - Spider Monkey Optimization (SMO)  

---

## ⚙️ Workflow

1. **Data Loading:** GEO dataset (GSE9476) is loaded and cleaned.  
2. **Hold-Out Split:** 70% training / 30% test.  
3. **Model 1 – SVM:** Using all genes (~12,625).  
4. **Model 2 – SVM + mRMR:** Select top 200 genes using mRMR.  
5. **Model 3 – RF + SMO Hybrid:** Random Forest ranks top 200 genes; SMO selects optimal subset (~87 genes).  
6. **Evaluation:** Accuracy, confusion matrices, and visualizations.  
7. **Visualizations:**  
   - Accuracy comparison bar chart  
   - Feature-selection funnel plot  
   - PCA scatterplot of SMO-selected genes  
   - Confusion matrices

---

## ✅ Results
- **SVM (All Genes):** ~78.95% accuracy  
- **SVM + mRMR (Top 200):** ~84.21% accuracy  
- **RF + SMO Hybrid:** ~94.74% accuracy with 87 genes selected  

---
