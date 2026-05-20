# =====================================================
# 1. LOAD LIBRARIES
# =====================================================

library(readxl)
library(survival)
library(survminer)
library(dplyr)


# =====================================================
# 2. IMPORT DATA
# =====================================================

data <- read_excel(
  "C:/Users/UTTATI/Downloads/Clinical IgaN.xlsx"
)


# =====================================================
# 3. SELECT VARIABLES
# =====================================================

cox_data <- data %>%
  
  select(
    
    `time(months)`,
    
    event,
    
    `Initial eGFR (mL/min/1.73m²)`,
    
    HTN,
    
    `Proteinuria 1`,
    
    `Proteinuria 2`,
    
    `Proteinuria 3`,
    
    `Proteinuria 4`
  )


# =====================================================
# 4. RENAME COLUMNS
# =====================================================

colnames(cox_data) <- c(
  
  "time",
  
  "event",
  
  "egfr",
  
  "HTN",
  
  "prot1",
  
  "prot2",
  
  "prot3",
  
  "prot4"
)


# =====================================================
# 5. CONVERT "NA" TEXT TO TRUE NA
# =====================================================

cox_data[cox_data == "NA"] <- NA


# =====================================================
# 6. REMOVE MISSING VALUES
# =====================================================

cox_complete <- na.omit(
  cox_data
)
dim(cox_complete)

# =====================================================
# 7. CONVERT VARIABLE TYPES
# =====================================================

cox_complete$time <- as.numeric(
  cox_complete$time
)

cox_complete$event <- as.numeric(
  cox_complete$event
)

cox_complete$egfr <- as.numeric(
  cox_complete$egfr
)

cox_complete$HTN <- factor(
  cox_complete$HTN
)

cox_complete$prot1 <- factor(
  cox_complete$prot1
)

cox_complete$prot2 <- factor(
  cox_complete$prot2
)

cox_complete$prot3 <- factor(
  cox_complete$prot3
)

cox_complete$prot4 <- factor(
  cox_complete$prot4
)


# =====================================================
# 8. CHECK DATA
# =====================================================

str(cox_complete)

summary(cox_complete)

colSums(is.na(cox_complete))


# =====================================================
# 9. FIT COX MODEL
# =====================================================

clinical_cox <- coxph(
  
  Surv(time, event) ~
    
    egfr +
    
    HTN +
    
    prot1 +
    
    prot2 +
    
    prot3 +
    
    prot4,
  
  data = cox_complete
)


# =====================================================
# 10. MODEL SUMMARY
# =====================================================

summary(clinical_cox)
AIC(clinical_cox)

# =====================================================
# 11. HAZARD RATIOS
# =====================================================

exp(
  coef(clinical_cox)
)


# =====================================================
# 12. CONFIDENCE INTERVALS
# =====================================================

exp(
  confint(clinical_cox)
)


# =====================================================
# 13. CONCORDANCE INDEX
# =====================================================

summary(
  clinical_cox
)$concordance


# =====================================================
# 14. ROYSTON-SAUERBREI R_D^2
# =====================================================

c_index <- summary(
  clinical_cox
)$concordance[1]

kappa <- sqrt(8/pi)

sigma2 <- pi^2 / 6

D <- qnorm(
  c_index
) * kappa

RD2 <- (D^2 / kappa^2) /
  
  (
    sigma2 +
      (D^2 / kappa^2)
  )

print(RD2)


# =====================================================
# 15. PROPORTIONAL HAZARDS TEST
# =====================================================

ph_test <- cox.zph(
  clinical_cox
)

print(ph_test)

plot(ph_test)


# =====================================================
# 16. FOREST PLOT
# =====================================================

ggforest(
  
  clinical_cox,
  
  data = cox_complete
)


# =====================================================
# 17. KM CURVE FOR HTN
# =====================================================

fit_htn <- survfit(
  
  Surv(time, event) ~ HTN,
  
  data = cox_complete
)

ggsurvplot(
  
  fit_htn,
  
  data = cox_complete,
  
  pval = TRUE,
  
  risk.table = TRUE,
  
  conf.int = FALSE
)


# =====================================================
# 18. KM CURVE FOR PROTEINURIA 1
# =====================================================

fit_prot1 <- survfit(
  
  Surv(time, event) ~ prot1,
  
  data = cox_complete
)

ggsurvplot(
  
  fit_prot1,
  
  data = cox_complete,
  
  pval = TRUE,
  
  risk.table = TRUE,
  
  conf.int = FALSE
)
AIC(clinical_cox)
