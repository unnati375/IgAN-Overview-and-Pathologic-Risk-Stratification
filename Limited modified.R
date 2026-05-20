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
  "C:/Users/UTTATI/Downloads/MESTC IgaN (1).xlsx"
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
    
    `Proteinuria 4`,
    
    M,
    
    E,
    
    S,
    
    T1,
    
    T2
  )


# =====================================================
# 4. RENAME VARIABLES
# =====================================================

colnames(cox_data) <- c(
  
  "time",
  
  "event",
  
  "egfr",
  
  "HTN",
  
  "prot1",
  
  "prot2",
  
  "prot3",
  
  "prot4",
  
  "M",
  
  "E",
  
  "S",
  
  "T1",
  
  "T2"
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
# 7. CONVERT VARIABLES TO CORRECT TYPES
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

cox_complete$M <- factor(
  cox_complete$M
)

cox_complete$E <- factor(
  cox_complete$E
)

cox_complete$S <- factor(
  cox_complete$S
)

cox_complete$T1 <- factor(
  cox_complete$T1
)

cox_complete$T2 <- factor(
  cox_complete$T2
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

mest_cox <- coxph(
  
  Surv(time, event) ~
    
    egfr +
    
    HTN +
    
    prot1 +
    
    prot2 +
    
    prot3 +
    
    prot4 +
    
    M +
    
    E +
    
    S +
    
    T1 +
    
    T2,
  
  data = cox_complete
)


# =====================================================
# 10. MODEL SUMMARY
# =====================================================

summary(
  mest_cox
)


# =====================================================
# 11. HAZARD RATIOS
# =====================================================

exp(
  coef(mest_cox)
)


# =====================================================
# 12. 95% CONFIDENCE INTERVALS
# =====================================================

exp(
  confint(mest_cox)
)


# =====================================================
# 13. AIC
# =====================================================

AIC(
  mest_cox
)


# =====================================================
# 14. CONCORDANCE INDEX
# =====================================================

summary(
  mest_cox
)$concordance


# =====================================================
# 15. Royston-Sauerbrei R_D²
# =====================================================

c_index <- summary(
  mest_cox
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
# 16. PROPORTIONAL HAZARDS TEST
# =====================================================

ph_test <- cox.zph(
  mest_cox
)

print(ph_test)

plot(ph_test)


# =====================================================
# 17. FOREST PLOT
# =====================================================

ggforest(
  
  mest_cox,
  
  data = cox_complete
)


# =====================================================
# 18. KAPLAN-MEIER CURVE FOR T1
# =====================================================

fit_T1 <- survfit(
  
  Surv(time, event) ~ T1,
  
  data = cox_complete
)

ggsurvplot(
  
  fit_T1,
  
  data = cox_complete,
  
  pval = TRUE,
  
  risk.table = TRUE,
  
  conf.int = FALSE,
  
  title = "Kaplan-Meier Curve by T1",
  
  xlab = "Follow-up Time",
  
  ylab = "Survival Probability"
)


# =====================================================
# 19. KAPLAN-MEIER CURVE FOR C2
# =====================================================

fit_C2 <- survfit(
  
  Surv(time, event) ~ C2,
  
  data = cox_complete
)

ggsurvplot(
  
  fit_C2,
  
  data = cox_complete,
  
  pval = TRUE,
  
  risk.table = TRUE,
  
  conf.int = FALSE,
  
  title = "Kaplan-Meier Curve by C2",
  
  xlab = "Follow-up Time",
  
  ylab = "Survival Probability"
)