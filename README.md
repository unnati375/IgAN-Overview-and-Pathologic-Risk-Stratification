# **IgA Nephropathy (IgAN) Risk Prediction in an Indian Cohort**

---

# **Abstract**

IgA Nephropathy (IgAN) is a chronic kidney disease caused by the accumulation of Immunoglobulin A (IgA) deposits in the kidney’s filtering units, known as glomeruli. Over time, these deposits trigger inflammation and reduce the kidney’s ability to effectively filter waste from the blood. Common clinical manifestations include:

- Hematuria (blood in urine)
- Proteinuria (protein in urine)
- Hypertension
- Progressive kidney dysfunction

While some patients remain clinically stable for years, others gradually progress to chronic kidney disease (CKD) or kidney failure.

Although IgAN has been extensively studied in European and East Asian populations, large-scale studies focusing on Indian patient cohorts remain limited. As a result, the disease progression patterns and associated risk factors within the Indian population are not yet fully understood.

This project aims to analyze the progression risk of IgAN in an Indian cohort from a medical institution in Kolkata using both clinical and histopathological variables. Statistical survival analysis techniques, particularly the **Cox Proportional Hazards Model**, are used to identify important predictors of disease progression. Additionally, feature attribution methods are applied to quantify the contribution of individual clinical and pathological factors toward patient risk prediction.

---

# **Problem Statement**

IgA Nephropathy (IgAN) is a progressive renal disease with highly variable clinical outcomes among patients. Although several internationally developed risk-prediction models exist, their applicability to Indian patient populations remains uncertain due to demographic, clinical, and pathological differences.

## **Key Challenges**

- Limited large-scale cohort studies on IgAN in India
- Differences in disease presentation and progression patterns across populations
- Lack of external validation of existing prediction models for South Asian cohorts
- Difficulty in accurately identifying high-risk patients at an early stage

This study addresses these challenges by developing and evaluating a survival-analysis-based framework tailored to an Indian IgAN cohort. The project combines:

- Clinical variables
- Laboratory findings
- MEST-C histopathological scores

to better understand disease progression patterns.

## **Objectives**

1. Identify significant clinical and histopathological predictors associated with IgAN progression.
2. Apply the Cox Proportional Hazards Model to estimate progression risk over time.
3. Evaluate whether incorporating MEST-C pathology scores improves predictive performance beyond clinical variables alone.
4. Quantify the contribution of individual risk factors using feature attribution techniques.

---

# **Research Questions**

This project investigates the following research questions:

- How effectively do clinical variables such as **eGFR**, **proteinuria**, and **hypertension** predict IgAN progression in an Indian patient cohort?

- Does the inclusion of **MEST-C histopathological scores** improve risk prediction compared to using clinical variables alone?

- Which individual features — including eGFR, hypertension status, proteinuria categories, and pathology scores — contribute most significantly to a patient’s predicted progression risk?

- Can survival-analysis-based prediction models provide clinically meaningful insights for early identification of high-risk IgAN patients in India?

---

# **Methodology Overview**

The study follows a survival analysis framework consisting of:

- Data preprocessing and missing value handling
- Exploratory Data Analysis (EDA)
- Clinical and pathological feature engineering
- Cox Proportional Hazards modeling
- Model comparison using AIC and concordance metrics
- Feature attribution and risk interpretation
- Survival probability estimation and visualization

---

# **Technologies & Tools**

- **R Programming**
- **Survival Analysis**
- **Cox Proportional Hazards Model**
- **MEST-C Scoring System**
- **Feature Attribution Techniques**
- **Statistical Data Visualization**

---

# **Expected Outcomes**

The project aims to:

- Identify major predictors of IgAN progression in Indian patients
- Evaluate the effectiveness of pathology-informed survival models
- Improve interpretability of clinical risk prediction
- Contribute toward population-specific evidence for kidney disease prognosis in India
