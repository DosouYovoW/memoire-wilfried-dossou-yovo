# R-squared
rsq <- function(y, y_hat) {
  sum((y_hat - mean(y))^2) / sum((y - mean(y))^2)
}


denman.beavers <- function(mat, maxit = 50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1, nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5 * (y + solve(z))
    z <- 0.5 * (z + solve(y))
    y <- y.temp
  }
  return(list(sqrt = y, sqrt.inv = z))
}


# Plot fixed effects and their confidence intervals and 
# importance of variables using t-values
nlmePlotConfTable <- function(mm, varCat, varCatOrder, paramCat,
                              conf_level = 0.95, 
                              pval.breaks = c(0, 0.05, 0.10, 1.00),
                              numApart = FALSE,
                              showRanef = FALSE,
                              varNamesFrom = NULL, varNamesTo = NULL,
                              varCatNamesTo = NULL, limits_df = NULL,
                              removeIntercept = FALSE,
                              varImpType = "t", # "t" or "cat", t is t-value, cat is correlation-adjusted t-scores
                              greyEnd = 0.7,
                              modify.param.names = NULL) {
  
  # Base table
  confTable <- as.data.frame(intervals(mm, level = conf_level)$fixed)
  
  if(showRanef) {
    for (i in seq_along(ranef(mm))) {
      confint_ranef <- t.test(ranef(mm)[[i]])$conf.int
      ranef_df <- data.frame(lower = confint_ranef[1], 
                             est. = mean(ranef(mm)[[i]]), 
                             upper = confint_ranef[2])
      rownames(ranef_df) <- paste0(names(ranef(mm))[[i]], '.ranef')
      confTable <- rbind(confTable, ranef_df)
    }
  }
  
  # p-values
  if(!showRanef) confTable$p <- summary(mm)$tTable[, 5]
  if(showRanef) confTable$p <- c(summary(mm)$tTable[, 5], rep(0, length(ranef(mm))))
  confTable$pCat <- cut(confTable$p, breaks = pval.breaks)
  confTable$pCat <- factor(sub(',', ' - ', confTable$pCat)) # prettier format
  confTable$test <- as.factor(ifelse(confTable[, 1] <= 0 & confTable[, 3] >= 0, "not significant", "significant"))
  
  # t-values without (t) or with (cat) adjustment for correlation
  if (varImpType == "t") {
    if(!showRanef) confTable$t <- summary(mm)$tTable[, 4]
    if(showRanef) confTable$t <- c(summary(mm)$tTable[, 4], rep(0, length(ranef(mm))))
    varImpXlab <- "t-value"
  } else if (varImpType == "cat") {
    # http://arxiv.org/pdf/0903.2003v4.pdf, eq 3.1, p.5
    # denman.beavers, custom function (defined) in this file to compute the sqrt of a matrix
    tcat <- denman.beavers(cov2cor(vcov(mm)))$sqrt.inv %*% summary(mm)$tTable[, 4]
    if(!showRanef) confTable$t <- tcat
    if(showRanef) confTable$t <- c(tcat, rep(0, length(ranef(mm))))
    varImpXlab <- "Correlation-adjusted t-scores"
  }
  
  # column for Variable category
  confTable$VarLong <- rownames(confTable) # to keep track of original names
  confTable$VarCat <- "Numeric"
  for (i in 1:length(varCat)) {
    confTable$VarCat[grepl(varCat[i], confTable$VarLong)] <- varCatNamesTo[i]
  }
  rownames(confTable) <- 1:nrow(confTable)
  confTable$VarCat <- as.factor(confTable$VarCat)
  confTable$VarCat <- factor(confTable$VarCat, levels = levels(confTable$VarCat)[varCatOrder])
  
  # Variable names
  confTable$Variable <- confTable$VarLong
  
  ## translate names
  if(!is.null(varNamesFrom)) {
    for (i in 1:length(varNamesFrom)) {
      confTable$Variable <- sub(pattern = varNamesFrom[i], 
                                replacement = varNamesTo[i], 
                                x = confTable$Variable)
    }
  }
  
  ## remove parameter names from variable names
  for (i in 1:length(paramCat)) {
    confTable$Variable <- sub(paste0(paramCat[i], '.'), '', confTable$Variable)
  }
  confTable$Variable <- factor(confTable$Variable)
  
  # column for Parameter category
  confTable$ParamCat <- NA
  for (i in 1:length(paramCat)) {
    filterMPC <- grepl(paramCat[i], confTable$VarLong)
    confTable$ParamCat[filterMPC] <- paramCat[i]
  }
  confTable$ParamCat <- as.factor(confTable$ParamCat)
  
  # limits 
  if(!is.null(limits_df)) {
    if(!is.na(limits_df[1, 2])) confTable$lower[confTable$ParamCat == paramCat[1] & confTable$lower < limits_df[1, 2]] <- NA
    if(!is.na(limits_df[1, 3])) confTable$upper[confTable$ParamCat == paramCat[1] & confTable$upper > limits_df[1, 3]] <- NA
    if(!is.na(limits_df[2, 2])) confTable$lower[confTable$ParamCat == paramCat[2] & confTable$lower < limits_df[2, 2]] <- NA
    if(!is.na(limits_df[2, 3])) confTable$upper[confTable$ParamCat == paramCat[2] & confTable$upper > limits_df[2, 3]] <- NA
    if(!is.na(limits_df[3, 2])) confTable$lower[confTable$ParamCat == paramCat[3] & confTable$lower < limits_df[3, 2]] <- NA
    if(!is.na(limits_df[3, 3])) confTable$upper[confTable$ParamCat == paramCat[2] & confTable$upper > limits_df[3, 3]] <- NA
  }
  
  ## change parameter names if specified
  if (!is.null(modify.param.names)) {
    levels(confTable$VarCat) <- modify.param.names
  }
  
  if(removeIntercept) {
    #confTable <- confTable[confTable$VarCat != 'Intercept', ]
    confTable <- confTable[confTable$Variable != '(Intercept)', ] # temporary hack
  }
  
  if(numApart) {
    fixefPar <- list()
    numPar <- droplevels(confTable[confTable$VarCat %in% c("Numeric"), ]) # "Intercept", 
    facPar <- droplevels(confTable[!(confTable$VarCat %in% c("Numeric")), ])
    levels(facPar$pCat) <- levels(confTable$pCat) # to assure that all p-values levels are in the legend
    levels(numPar$pCat) <- levels(confTable$pCat) # to assure that all p-values levels are in the legend
    fixefPar[[1]] <- ggplot(numPar, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y') +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 0.6) +
      geom_point(aes(colour=pCat), size = 1, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value', drop = FALSE) +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
    fixefPar[[2]] <- ggplot(facPar, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y') +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 0.6) +
      geom_point(aes(colour=pCat), size = 1, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value', drop = FALSE) +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
  } else {
    fixefPar <- ggplot(confTable, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 0.6) +
      geom_point(aes(colour=pCat), size = 1, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value') +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
    
    if(length(unique(confTable$ParamCat)) == 1) {
      fixefPar <- fixefPar + facet_grid(VarCat ~ ., scale = 'free', space = 'free_y')
    } else {
      fixefPar <- fixefPar + facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y')
    }
    
  }
  
  varImp <- ggplot(confTable, aes(x = abs(t), y = Variable)) +
    geom_segment(aes(x = 0, xend=abs(t), yend = Variable,
                     colour=pCat), lwd = 2) +
    scale_colour_grey(start = 0, end = greyEnd, name = 'p-value') +
    xlab(varImpXlab) +
    ylab('') +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0))
  
  if(length(unique(confTable$VarCat)) == 1) {
    varImp <- varImp + facet_grid(. ~ ParamCat, scale = 'free', space = 'free_y')
  } else {
    varImp <- varImp + facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y')
  }
  
  return(list(fixefPar, varImp))
}


prepare_nlme = function(X, variable) {
  
  ## scale
  varnumcols = sapply(X, class) %in% c("numeric", "integer") & 
    colnames(X) %in% variable
  X[, varnumcols] = apply(X[, varnumcols], 2, scale)
  
  ## Valeurs de départ
  start_list = list()
  for (i in 1:length(variable)) {
    if (is.factor(X[variable[i]][[1]])) {
      start_list[[i]] = rep(0, length(levels(X[variable[i]][[1]]))-1)
    } else {
      start_list[[i]] = 0
    }
  }
  start_vector = unlist(start_list)
  
  ## Variables explicatives
  rhs = paste(variable, collapse= " + ")
  
  return(list(start = start_vector, rhs = rhs, scaled = X))
}


# Yield response Mitscherlich
mitschFunc <- function(dose, asymptote, taux, environnement, rateExp = 'exp')  {
  if (rateExp == 'exp') {
    asymptote * (1 - exp(-exp(taux) * (dose + environnement)))
  } else {
    asymptote * (1 - exp(-taux * (dose + environnement)))
  }
} 

## Mitscherlich prediction
mitschPred <- function(mm, newdata, rhs, col_dose, ranEf = 0) {
  library(stringi)
  # parameters
  ## collect parameters from the model
  parameter <- summary(mm)$tTable
  
  ## collect specific mean and standard error parameters for Asymptote, Taux and Environnement
  asymParam <- parameter[stri_detect_fixed(rownames(parameter), "Asym"), 1]
  tauxParam <- parameter[stri_detect_fixed(rownames(parameter), "Rate"), 1] 
  enviParam <- parameter[stri_detect_fixed(rownames(parameter), "Envi"), 1]
  
  ## the linear combination of site specific condition (model.matrix) used to compute Mitscherlich parameters
  modmat <- model.matrix(as.formula(paste('~', paste(rhs, collapse = '+'))), data = newdata)
  asymMitsch <- modmat %*% asymParam
  tauxMitsch <- modmat %*% tauxParam
  enviMitsch <- modmat %*% enviParam
  
  # predict
  doseMitsch <- newdata[[col_dose]]
  
  # Build formula according to specified constrains
  asymTag <- asymMitsch[1] + ranEf
  envTag <- enviMitsch[1]
  rateTag <- tauxMitsch[1]
  
  pred <- asymTag * (1 - exp(-rateTag * (doseMitsch + envTag)))
  
  # output
  return(list(pred = pred, 
              metaParam = list(A = asymMitsch[1], R = tauxMitsch[1], E = enviMitsch[1])))
}

# Optimal dose Mitscherlich
doseOptMitsherlich <- function(asymptote, taux, environnement, 
                               prix_dose, prix_vendable, interets = 0, periode = 0, 
                               rateExp = 'exp') {
  if (rateExp == 'exp') {
    -exp(-taux) * (log(prix_dose * (1 + interets)**periode / (asymptote * prix_vendable)) + exp(taux) * environnement - taux)
  } else {
    (log((asymptote * taux * prix_vendable) / (prix_dose * (interets + 1)**periode)) - environnement * taux) / taux 
  }
}


