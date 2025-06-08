#############
###  function to get the overall meta-analytic proportion and the effect of Time (in 4  periods)
############

# x = character string of the variable
# Ntot = character string of the variable with total sample size
# name = character string of the name that will appear in the table
# time = character string of the period variable

PROP1 <- function(x, Ntot, name, time)
{
  
  dat <- data.frame(count=as.integer(DATA[,x]),
                    Ntot=DATA[,Ntot], 
                    t=DATA[,time]) %>%
    mutate(prop = count/Ntot,
           wt = log(Ntot))
  dat <- na.omit(dat)
  
  N_tot   <-  nrow(dat)
  N_time <- table(dat$t)
  
  #global model 
  overall_mod <- glm(prop ~ 1, family = "binomial", weights = wt, data = dat)
  ov_prop <- as.data.frame(emmeans(overall_mod, ~ 1, type = "response")) %>%
    mutate_if(is.numeric, function(x){round(x*100, 0)}) %>% 
    mutate(prop_CI = paste0(prob, " (", asymp.LCL,"-", asymp.UCL,")"))
  
  #model with time
  time_mod <- glm(prop ~ t, family = "binomial", weights = wt, data = dat)
  time_prop <- as.data.frame(emmeans(time_mod, ~ t, type = "response")) %>%
    mutate_if(is.numeric, function(x){round(x*100, 0)}) %>% 
    mutate(prop_CI = paste0(prob, " (", asymp.LCL,"-", asymp.UCL,")"))
  
  #significance
  pval<- ifelse(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2] < 0.0001, "<0.0001", round(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2],4))
  
  res <- c(name, N_tot, ov_prop$prop_CI, 
           N_time[1], time_prop$prop_CI[1],
           N_time[2], time_prop$prop_CI[2],
           N_time[3], time_prop$prop_CI[3],
           N_time[4], time_prop$prop_CI[4], pval)
  
  return(res)
}

#########
### function to calculate proporiton of studies within each time period for continuous variables that were transformed in a cutoff
#########

# x = character string of the variable
# name = character string of the name that will appear in the table
# time = character string of the time period variable
# method = either "chisq" or "fisher" depending on the type of the test wanted

CAT<-function(x, name, time, method)
{
      tabtot <- table(DATA[,x])
      sumtot <- sum(tabtot)
      tabtime <- table(DATA[,x], DATA[,time])
      sumtime <- apply(tabtime, 2, sum)
  
      proptot <- tabtot[2]/sum(tabtot)
      CItot <- exactci(tabtot[2], sum(tabtot), conf.level = 0.95)
      
      proptime <- tabtime[2,]/sumtime
      CItime <- map2(tabtime[2,], sumtime,
                    function(x, y){exactci(x, y, conf.level = 0.95)})
      propCItime <- map2(proptime, CItime, 
                    function(x,y){paste0(round(x*100,0)," (", 
                                        round(y[["conf.int"]][1]*100,0),"-",
                                        round(y[["conf.int"]][2]*100,0),")")})

  if(method == "fisher"){
        p <- fisher.test(tabtime)
  } else if(method == "chisq"){
        p<-chisq.test(tabtime)
  }
      
      N_propCItime <- NULL
      for(i in 1:length(propCItime)){
        N_propCItime <- c(N_propCItime, c(sumtime[i], propCItime[[i]]))
      }
      
      res=c(name, 
            sumtot,
            paste0(round(proptot*100,0)," (",round(CItot$conf.int[1]*100,0),"-", round(CItot$conf.int[2]*100,0),")"),
            N_propCItime,
            round(p$p.value[1],4))
      return(res)
}

##########
### function to calculate overall effect size for death or readmission rates (output has results + model to do model checks)
##########

# var = character string of the variable to calculate the effect size for

rate_ES <- function(var){
        dat <- na.omit(DATA[,c("number_follow_up", var)])
        colnames(dat)<-c("Ni","Ei")
  
        dat_es <- escalc(
          xi= Ei,
          ni =Ni,
          data = dat, 
          measure = "PLO",
          to="if0all")

  
  #global effect
          ov_mod <- rma(yi=yi, vi=vi, data=dat_es, method = "REML")
  
  #prediction
          pred_overall <- predict(ov_mod,  transf=transf.ilogit, addx=T)
          pred_overall <- as.data.frame(pred_overall)
          
  #extraction
          Eov  <- round(pred_overall$pred*100, 0)
          CLov <- round(pred_overall$ci.lb*100, 0)
          CUov <- round(pred_overall$ci.ub*100, 0)
          I2ov <- round(ov_mod$I2, 2)
          
          summ_rate <- data.frame(N = nrow(dat),
                                  "Estimate (%)" = Eov,
                                  "95% CI" = paste0("(",CLov," - ", CUov, ")"),
                                  "I^2" =  I2ov, check.names = FALSE)
  
  return(list(summ_rate, ov_mod))
}


###########
### function to get weighted logistic regression with continuous X variable and proportion Y variable and ggplot showing the predicted relationship
###########

# x = X variable as character
# y = Y variable as character
# Xpred = vector of X values to get predictions for
# labX = X label as character
# labY = Y label as character
# limY = axis limits for Y axis, vector of 2 values
# title = plot title as character
# ypos = y relative position of text with OR and p-value


wtlogis_plot <- function(NB, x, y, Xpred, labX, labY, limY, title, ypos, unit){
  
  dat <- data.frame(count=as.integer(DATA[,y]),
                    Ntot=DATA[,NB], 
                    t=DATA[,x]) %>%
    mutate(prop = count/Ntot,
           perc = prop*100,
           wt = log(Ntot),
           t_div = t/10)
  dat <- na.omit(dat)
  
  
  #global model 
  overall_mod <- glm(prop ~ 1, family = "binomial", weights = wt, data = dat)
  
  #model with time
  time_mod <- glm(prop ~ t, family = "binomial", weights = wt, data = dat)
  time_rg <- ref_grid(time_mod, at = list(t = Xpred))
  
  time_prop <- as.data.frame(emmeans(time_rg, ~ t, type = "response")) %>%
    mutate_at(c(2,5:6), function(x){x*100})
  begend <- round(c(time_prop$prob[1], last(time_prop$prob)),0)
  
  # model with time divided by 10 to get OR for 10yr increment
  timediv_mod <- glm(prop ~ t_div, family = "binomial", weights = wt, data = dat)
  ORdiv <- tidy(timediv_mod, conf.int = TRUE) %>%
    mutate(OR_CI = paste0(round(exp(estimate),2), " (", 
                          round(exp(conf.low),2),"-", 
                          round(exp(conf.high),2), ")"))
  
  #significance
  pval<- ifelse(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2] < 0.001, 
                "<0.001", 
                round(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2],3))
  
  OR_fig <- c(ORdiv$OR_CI[2], pval)
  
  #sample size specific data points
  dat$size <-  2 + 20* ((1 / sqrt(dat$Ntot)) - min((1 / sqrt(dat$Ntot)),na.rm=T))/(max((1 / sqrt(dat$Ntot)),na.rm=T) - min((1 / sqrt(dat$Ntot)),na.rm=T))
  
  #Graph
  mods_plot <- ggplot(data=dat, aes(x= t, y= perc)) +
    geom_line(data=time_prop, aes(x= t, y= prob), size=2, colour="blue") +
    geom_line(data=time_prop, aes(x= t, y= asymp.LCL), size=1,linetype = 2, colour="blue") +
    geom_line(data=time_prop, aes(x= t, y= asymp.UCL), size=1,linetype = 2, colour="blue") +
    geom_jitter(size = dat$wt, shape=1, stroke=1,width = 0.3) +
    scale_color_manual(values=c("blue")) +
    guides(size=FALSE) + ylab(labY) + xlab(labX) +
    theme_classic() +
    theme(axis.title.x = element_text( face="bold", size=30),
          axis.title.y = element_text( face="bold", size=30),
          axis.text.x = element_text(face="bold", size=30),
          axis.text.y = element_text(face="bold", size=30),
          plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
          legend.position = "none") +
    ylim(limY) + ggtitle(title)
  
  OR_p_plot <- paste0("OR (CI):", OR_fig[1],"\n", "N=", nrow(dat),"; p=", OR_fig[2],"\nfrom ", begend[1], " to ", begend[2], " ", unit)
  grob_ORp <- grobTree(textGrob(OR_p_plot, x=0.02,  y= ypos, hjust=0,
                                gp=gpar(col="black", fontsize=25, fontface="bold")))
  
  mods_plot <- mods_plot + annotation_custom(grob_ORp)
  
  return(list(time_mod, mods_plot, OR_fig, timediv_mod))
}

###########
### function to get weighted linear regression with continuous X variable and continuous Y variable and ggplot showing the predicted relationship
###########

# x = X variable as character
# y = Y variable as character
# Xpred = vector of X values to get predictions for
# labX = X label as character
# labY = Y label as character
# limY = axis limits for Y axis, vector of 2 values
# title = plot title as character
# ypos = y relative position of text with OR and p-value


wtlin_plot <- function(NB, x, y, Xpred, labX, labY, limY, title, ypos, unit){
  
  dat <- data.frame(y=as.integer(DATA[,y]),
                    Ntot=DATA[,NB], 
                    t=DATA[,x]) %>%
    mutate(wt = log(Ntot),
           t_div = t/10)
  dat <- na.omit(dat)
  
  
  #global model 
  overall_mod <- lm(y ~ 1, weights = wt, data = dat)
  
  #model with time
  time_mod <- lm(y ~ t, weights = wt, data = dat)
  time_rg <- ref_grid(time_mod, at = list(t = Xpred))
  
  time_pred <- as.data.frame(emmeans(time_rg, ~ t))
  begend <- round(c(time_pred$emmean[1], last(time_pred$emmean)),0)
  
  # model with time divided by 10 to get beta for 10yr increment
  timediv_mod <- lm(y ~ t_div, weights = wt, data = dat)
  betadiv <- tidy(timediv_mod, conf.int = TRUE) %>%
    mutate(beta_CI = paste0(round(estimate,2), " (", 
                            round(conf.low,2),"-", 
                            round(conf.high,2), ")"))
  
  #significance
  pval<- ifelse(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2] < 0.001, 
                "<0.001", 
                round(anova(overall_mod,time_mod, test= "Chisq")$`Pr(>Chi)`[2],3))
  
  beta_fig <- c(betadiv$beta_CI[2], pval)
  
  #sample size specific data points
  dat$size <-  2 + 20* ((1 / sqrt(dat$Ntot)) - min((1 / sqrt(dat$Ntot)),na.rm=T))/(max((1 / sqrt(dat$Ntot)),na.rm=T) - min((1 / sqrt(dat$Ntot)),na.rm=T))
  
  #Graph
  mods_plot <- ggplot(data=dat, aes(x= t, y= y)) +
    geom_line(data=time_pred, aes(x= t, y= emmean), size=2, colour="blue") +
    geom_line(data=time_pred, aes(x= t, y= lower.CL), size=1,linetype = 2, colour="blue") +
    geom_line(data=time_pred, aes(x= t, y= upper.CL), size=1,linetype = 2, colour="blue") +
    geom_jitter(size = dat$wt, shape=1, stroke=1,width = 0.3) +
    scale_color_manual(values=c("blue")) +
    guides(size=FALSE) + ylab(labY) + xlab(labX) +
    theme_classic() +
    theme(axis.title.x = element_text( face="bold", size=30),
          axis.title.y = element_text( face="bold", size=30),
          axis.text.x = element_text(face="bold", size=30),
          axis.text.y = element_text(face="bold", size=30),
          plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
          legend.position = "none") +
    ylim(limY) + ggtitle(title)
  
  beta_p_plot <- paste0("beta (CI):", beta_fig[1],"\n", "N=", nrow(dat),"; p=", beta_fig[2],"\nfrom ", begend[1], " to ", begend[2], " ", unit)
  grob_betap <- grobTree(textGrob(beta_p_plot, x=0.02,  y= ypos, hjust=0,
                                  gp=gpar(col="black", fontsize=25, fontface="bold")))
  
  mods_plot <- mods_plot + annotation_custom(grob_betap)
  
  return(list(time_mod, mods_plot, beta_fig, timediv_mod))
}



###########
### function to get OR and CI from the model with 10 year increment
##########

# obj = model with x variable divided by 10
# p = p-value from the LRT

OR_div <- function(obj, p){
        OR_CI <- paste0(round(exp(obj$beta[2]),2), " (", 
                        round(exp(obj$ci.lb[2]),2),"-", 
                        round(exp(obj$ci.ub[2]),2), ")")
        pval_clean <- ifelse(p[["pval"]] < 0.001, "p<0.001", 
                             paste0("p=",round(p[["pval"]],3)))
  
        OR_p <- c(OR_CI, p, pval_clean)
  return(OR_p)
}


###########
### function to get meta-regression model with continuous X variable and proportion Y variable and ggplot showing the predicted relationship
###########

# x = X variable as character
# y = Y variable as character
# Xpred = vector of X values to get predictions for
# labX = X label as character
# labY = Y label as character
# limY = axis limits for Y axis, vector of 2 values
# title = plot title as character
# ypos = y relative position of text with OR and p-value


MA_prop_plot <- function(NB, x, y, Xpred, labX, labY, limY, title, ypos){
  
              prop_dat <- na.omit(DATA[,c(NB, y, x)])
              colnames(prop_dat)<-c("Ni","Ei","X")
              prop_dat$prop <- (prop_dat$Ei/prop_dat$Ni) * 100
              prop_dat$X_div <- prop_dat$X/10
              
              es_prop_dat <- escalc(
                xi= Ei,
                ni =Ni,
                data = prop_dat, 
                measure = "PLO",
                to="if0all")
  
  #global effect
              mod_overall <- rma(yi=yi, vi=vi, data=es_prop_dat, method="REML")
  
  #moderator
              mod_mods <- rma(yi=yi,mods = ~X, vi=vi, data=es_prop_dat, method = "REML")
              mod_mods_div <- rma(yi=yi,mods = ~X_div, vi=vi, data=es_prop_dat, method = "REML")

  #significance test
              mod_overall_ML <- rma(yi=yi, vi=vi, data=es_prop_dat, method = "ML") 
              mod_mods_ML <- rma(yi=yi,mods = ~X, vi=vi, data=es_prop_dat, method = "ML")
              pval <- anova(mod_overall_ML, mod_mods_ML)[6]
              
              OR_fig <- OR_div(mod_mods_div, pval)
  
  #Prediction
              pred_mods <- as.data.frame(predict(mod_mods, newmods = Xpred, transf=transf.ilogit, addx=T))
              begend <- round(c(pred_mods$pred[1], last(pred_mods$pred))*100,0)
  
  #sample size specific data points
              es_prop_dat$size <-  2 + 20* ((1 / sqrt(es_prop_dat$vi)) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))/(max((1 / sqrt(es_prop_dat$vi)),na.rm=T) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))
  
  #Graph
            mods_plot <- ggplot(data=es_prop_dat, aes(x= X, y= prop)) +
              geom_line(data=pred_mods, aes(x= Xpred, y= pred*100), size=2, colour="blue") +
              geom_line(data=pred_mods, aes(x= Xpred, y= ci.lb*100), size=1,linetype = 2, colour="blue") +
              geom_line(data=pred_mods, aes(x= Xpred, y= ci.ub*100), size=1,linetype = 2, colour="blue") +
              geom_jitter(size = es_prop_dat$size, shape=1, stroke=1,width = 0.3) +
              scale_color_manual(values=c("blue")) +
              guides(size=FALSE) + ylab(labY) + xlab(labX) +
              theme_classic() +
              theme(axis.title.x = element_text( face="bold", size=30),
                    axis.title.y = element_text( face="bold", size=30),
                    axis.text.x = element_text(face="bold", size=30),
                    axis.text.y = element_text(face="bold", size=30),
                    plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
                    legend.position = "none") +
              ylim(limY) + ggtitle(title)
  
              OR_p_plot <- paste0("OR (CI):", OR_fig[1],"\n", "N=", nrow(es_prop_dat),"; ", OR_fig[3],"\nfrom ", begend[1]," to ", begend[2],"%")
              grob_ORp <- grobTree(textGrob(OR_p_plot, x=0.02,  y= ypos, hjust=0,
                                            gp=gpar(col="black", fontsize=25, fontface="bold")))
              
              mods_plot <- mods_plot + annotation_custom(grob_ORp)
  
  return(list(mod_mods, mods_plot, OR_fig, mod_mods_div))
}

#### very similar function but non-linear trends without adjustment on continent

MA_prop_nonlin_plot <- function(NB, x, y, Xpred, labX, labY, limY, title, ypos, knots){
  
  prop_dat <- na.omit(DATA[,c(NB, y, x)])
  colnames(prop_dat)<-c("Ni","Ei","X")
  prop_dat$prop <- (prop_dat$Ei/prop_dat$Ni) * 100
  prop_dat$X_div <- prop_dat$X/10
  
  es_prop_dat <- escalc(
    xi= Ei,
    ni =Ni,
    data = prop_dat, 
    measure = "PLO",
    to="if0all")
  
  #global effect
  mod_overall <- rma(yi=yi, vi=vi, data=es_prop_dat, method="REML")
  
  #moderator non-linear
  mod_mods <- rma(yi=yi,mods = ~X + rcspline.eval(X, nk = knots, inclx = FALSE), 
                      vi=vi, data=es_prop_dat, method = "REML")
  mod_mods_div <- rma(yi=yi,mods = ~X + rcspline.eval(X_div, nk = knots, inclx = FALSE), 
                      vi=vi, data=es_prop_dat, method = "REML")
  
  #significance test for non-linear component
  mod_lin_ML <- rma(yi=yi, vi=vi, mods= ~X, data=es_prop_dat, method = "ML") 
  mod_nonlin_ML <- rma(yi=yi,mods = ~X + rcspline.eval(X, nk = knots,inclx = FALSE), 
                       vi=vi, data=es_prop_dat, method = "ML")
  pval <- anova(mod_lin_ML, mod_nonlin_ML)[6]
  
  OR_fig <- OR_div(mod_mods_div, pval)
  
  knots_pos <- attr(rcspline.eval(es_prop_dat$X, nk = 3, inclx = TRUE), "knots")
  
  #Prediction
  mod_pred <- rma(yi=yi,mods = ~rcspline.eval(X, nk = knots, inclx = TRUE), vi=vi, data=es_prop_dat, method = "REML")
  pred_mods <- as.data.frame(predict(mod_pred, newmods = cbind(rcspline.eval(Xpred, knots_pos, inclx = TRUE)),
                                     transf=transf.ilogit, addx=T))
  
  #sample size specific data points
  es_prop_dat$size <-  2 + 20* ((1 / sqrt(es_prop_dat$vi)) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))/(max((1 / sqrt(es_prop_dat$vi)),na.rm=T) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))
  
  #Graph
  mods_plot <- ggplot(data=es_prop_dat, aes(x= X, y= prop)) +
    geom_line(data=pred_mods, aes(x= Xpred, y= pred*100), size=2, colour="blue") +
    geom_line(data=pred_mods, aes(x= Xpred, y= ci.lb*100), size=1,linetype = 2, colour="blue") +
    geom_line(data=pred_mods, aes(x= Xpred, y= ci.ub*100), size=1,linetype = 2, colour="blue") +
    geom_jitter(size = es_prop_dat$size, shape=1, stroke=1,width = 0.3) +
    scale_color_manual(values=c("blue")) +
    guides(size=FALSE) + ylab(labY) + xlab(labX) +
    theme_classic() +
    theme(axis.title.x = element_text( face="bold", size=30),
          axis.title.y = element_text( face="bold", size=30),
          axis.text.x = element_text(face="bold", size=30),
          axis.text.y = element_text(face="bold", size=30),
          plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
          legend.position = "none") +
    ylim(limY) + ggtitle(title)
  
  #OR_p_plot <- paste0("OR (CI):", OR_fig[1],"\n", "N=", nrow(es_prop_dat),"; ", OR_fig[3])
  # grob_ORp <- grobTree(textGrob(OR_p_plot, x=0.02,  y= ypos, hjust=0,
  #                             gp=gpar(col="black", fontsize=30, fontface="bold")))
  
  grob_p <- grobTree(textGrob(paste0(" p non-linear component= ", round(pval$pval,3)), x=0.02,  y= ypos, hjust=0,
                              gp=gpar(col="black", fontsize=20, fontface="bold")))
  
  mods_plot <- mods_plot + annotation_custom(grob_p)
  
  return(list(mod_mods, mods_plot, mod_mods_div, pred_mods))
}




#########
### function to get OR for subgroup variables with 2 categories
#########

# x = character string for the x variable
# y = character string for the y variable
# group = chracter string for the group variable
# newref = character string for the level of the group variable that is not the reference level

OR_group_fun <- function(x, y, group, newref){
  
          OR_data <- data.frame(Group = rep(NA,2), N = rep(NA,2), OR = rep(NA,2), lowCI = rep(NA,2), upCI = rep(NA,2), p = rep(NA,2), I2 = rep(NA,2))
  
          group_dat <- na.omit(DATA[,c("number_follow_up", y, x, group)])
          colnames(group_dat)<-c("Ni","Ei","X","group")
          group_dat$X_div <- group_dat$X/10
          
          OR_data$Group <- levels(group_dat$group)
          OR_data$N <- as.integer(table(group_dat$group))
          
          es_group_dat <- escalc(
            xi= Ei,
            ni =Ni,
            data = group_dat, 
            measure = "PLO",
            to="if0all")
  
          inter_mod1 <- rma(yi=yi,mods = ~ X_div * group, vi=vi, data=es_group_dat, method = "REML")
          inter_mod2 <- rma(yi=yi,mods = ~ X_div * I(relevel(group, newref)) , vi=vi, data=es_group_dat, method = "REML")
          
          #significance test
          inter_mod_full <- rma(yi=yi, mods =  ~ X_div * group, vi=vi, data=es_group_dat, method = "ML") 
          inter_mod_red <- rma(yi=yi, mods =  ~ X_div + group, vi=vi, data=es_group_dat, method = "ML")
          pval <- anova(inter_mod_full, inter_mod_red)[6]
          OR_data$p[2] <- ifelse(pval$pval < 0.001, "<0.001", round(pval$pval,3))
          
          OR_data$OR <- c(exp(inter_mod1$b[2]), exp(inter_mod2$b[2]))
          OR_data$lowCI <- c(exp(inter_mod1$ci.lb[2]), exp(inter_mod2$ci.lb[2]))
          OR_data$upCI <- c(exp(inter_mod1$ci.ub[2]), exp(inter_mod2$ci.ub[2]))
          OR_data$I2 <- round(inter_mod1$I2,2)
  
  return(OR_data)
}


########
### function to get meta-regression model with continuous X1 and X2 variables and 2 group variables and proportion Y variable and ggplot showing the predicted relationship with x1 for a median value of x2
########

# x1 = X1 variable as character
# x2 = X2 variable as character
# group1 = group1 variable (2-level factor) as character
# group2 = group2 variable (2-level factor) as character
# y = Y variable as character
# x1pred = vector of values for which to predict x1
# x2val = median value of X2 for which to get predictions
# labX = X label as character
# labY = Y label as character
# limY = axis limits for Y axis, vector of 2 values
# title = plot title as character
# ypos = y relative for text of OR and p-value

MA_prop_plot_adj <- function(NB, x1, x2, group1, group2, y, x1pred, x2val, labX, labY, limY, title, ypos){
  
          prop_dat <- na.omit(DATA[,c(NB, y, x1, x2, group1, group2)])
          colnames(prop_dat)<-c("Ni","Ei","X1", "X2", "group1", "group2")
          prop_dat$prop <- (prop_dat$Ei/prop_dat$Ni) * 100
          prop_dat$X1_div <- prop_dat$X1/10
          
          es_prop_dat <- escalc(
            xi= Ei,
            ni =Ni,
            data = prop_dat, 
            measure = "PLO",
            to="if0all")
  
  # adapts full model formula depending on number of studies in each group
    if(all(table(prop_dat$group1) > 10) & all(table(prop_dat$group2) > 10)){
            form <- ~X1 + X2 + group1 + group2 # add group1 & group2 when both have N>10
            form1 <- ~X1_div + X2 + group1 + group2
            form_red <- ~X2 + group1 + group2
    } else if(!all(table(prop_dat$group1) > 10) & !all(table(prop_dat$group2) > 10)){
            form <- ~X1 + X2
            form1 <- ~X1_div + X2
            form_red <- ~X2
    } else if(all(table(prop_dat$group1) > 10) & !all(table(prop_dat$group2) > 10)){
            form <- ~ X1 + X2 + group1
            form1 <- ~X1_div + X2 + group1
            form_red <- ~X2 + group1
    } else {form <-  ~ X1 + X2 + group2
            form1 <- ~ X1_div + X2 + group2
            form_red <- ~X2 + group2
            pred_mat <- as.matrix(expand.grid(X1 = x1pred, X2 = x2val, group2 = c(0,1)))
    }
  
            form_nogp <- ~ X1 + X2
            form_nogp1 <- ~ X1_div + X2
            form_nogp_red <- ~ X2
            pred_mat <- as.matrix(expand.grid(X1 = x1pred, X2 = x2val))
            
    #full model
            mod_full <- rma(yi=yi, mods = form, vi=vi, data=es_prop_dat, method="REML")
            VIF <- vif(mod_full)
            mod_full_div <- rma(yi=yi, mods = form1, vi=vi, data=es_prop_dat, method="REML")
            
    #reduced model
            mod_red <- rma(yi=yi, mods = form_red, vi=vi, data=es_prop_dat, method="REML")
  
    #significance test
            mod_full_ML <- rma(yi=yi, mods = form, vi=vi, data=es_prop_dat, method = "ML") 
            mod_red_ML <- rma(yi=yi, mods = form_red, vi=vi, data=es_prop_dat, method = "ML")
            pval <- anova(mod_full_ML, mod_red_ML)[6]
            pval_clean <- ifelse(pval$pval < 0.001, "p<0.001", paste0("p=",round(pval$pval,4)))
            
            OR_gp <- OR_div(mod_full_div, pval)
  
  
    ### models with no group (just Year as moderator)
      #full model
            mod_nogp_full <- rma(yi=yi, mods = form_nogp, vi=vi, data=es_prop_dat, method="REML")
            mod_nogp_full_div <- rma(yi=yi, mods = form_nogp1, vi=vi, data=es_prop_dat, method="REML")
  
      #reduced model
            mod_nogp_red <- rma(yi=yi, mods = form_nogp_red, vi=vi, data=es_prop_dat, method="REML")
  
    #significance test
            mod_nogp_full_ML <- rma(yi=yi, mods = form_nogp, vi=vi, data=es_prop_dat, method = "ML") 
            mod_nogp_red_ML <- rma(yi=yi, mods = form_nogp_red, vi=vi, data=es_prop_dat, method = "ML")
            pval_nogp <- anova(mod_nogp_full_ML, mod_nogp_red_ML)[6]
            pval_nogp_clean <- ifelse(pval_nogp$pval < 0.001, "p<0.001", paste0("p=",round(pval_nogp$pval,4)))
            
            OR_nogp <- OR_div(mod_nogp_full_div, pval_nogp)
  
    #Prediction
            pred_mods <- as.data.frame(predict(mod_nogp_full, newmods = pred_mat, transf=transf.ilogit, addx=T))
            pred_mods_avg <- pred_mods %>% 
              group_by(X.X1) %>% 
              summarise(mean_pred = mean(pred),
                        mean_ci.lb = mean(ci.lb),
                        mean_ci.ub = mean(ci.ub))
            begend <- round(c(pred_mods_avg$mean_pred[1], last(pred_mods_avg$mean_pred))*100,0)
  
    #sample size specific data points
              es_prop_dat$size <-  2 + 20* ((1 / sqrt(es_prop_dat$vi)) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))/(max((1 / sqrt(es_prop_dat$vi)),na.rm=T) - min((1 / sqrt(es_prop_dat$vi)),na.rm=T))
  
    #Graph
        mods_plot <- ggplot(data=es_prop_dat, aes(x= X1, y= prop)) +
              geom_line(data=pred_mods_avg, aes(x= X.X1, y= mean_pred*100), size=2, colour="blue") +
              geom_line(data=pred_mods_avg, aes(x= X.X1, y= mean_ci.lb*100), size=1,linetype = 2, colour="blue") +
              geom_line(data=pred_mods_avg, aes(x= X.X1, y= mean_ci.ub*100), size=1,linetype = 2, colour="blue") +
              geom_jitter(size = es_prop_dat$size, shape=1, stroke=1,width = 0.3) +
              scale_color_manual(values=c("blue")) +
              guides(size=FALSE) + ylab(labY) + xlab(labX) +
              theme_classic() +
              theme(axis.title.x = element_text( face="bold", size=30),
                    axis.title.y = element_text( face="bold", size=30),
                    axis.text.x = element_text(face="bold", size=30),
                    axis.text.y = element_text(face="bold", size=30),
                    plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
                    legend.position = "none") +
              ylim(limY) + ggtitle(title)
            
            OR_p_plot <- paste0("OR (CI):", OR_nogp[1],"\n","N=", nrow(es_prop_dat), "; ", OR_nogp[3], "\nfrom ", begend[1]," to ", begend[2],"%")
            
            grob <- grobTree(textGrob(OR_p_plot, x=0.02,  y=ypos, hjust=0,
                                      gp=gpar(col="black", fontsize=25, fontface="bold")))
            
            mods_plot <- mods_plot + annotation_custom(grob)
  
  return(list(mod_full_div, VIF, mods_plot, OR_gp, OR_nogp))
}


########
### function to get OR of time and p-values for subgroup variables with 2 categories 
########

# NB = character string for the variable with sample size
# t = character string for the time variable
# y = character string for the variable with the number of events
# group = character string for the group variable
# newref = character string for the level of the group variable that is not the reference level
# var = character string for the name of the variable in the table

sub_OR_group_fun <- function(NB, t, y, group, newref, var){
  
          OR_data <- data.frame(Variable = c(var, "", ""), Group = rep(NA,3), N = rep(NA,3), OR = rep(NA,3), lowCI = rep(NA,3), upCI = rep(NA,3), p = rep(NA, 3), "p inter" = rep(NA,3), check.names = FALSE)
          
          group_dat <- na.omit(DATA[,c(NB, y, t, group)])
          colnames(group_dat)<-c("Ni","Ei","time","group")
          group_dat$time_div <- group_dat$time/10
          
          OR_data$Group <- c("", levels(group_dat$group))
          OR_data$N <- c("", table(group_dat$group))
          
          es_group_dat <- escalc(
            xi= Ei,
            ni =Ni,
            data = group_dat, 
            measure = "PLO",
            to="if0all")
  
  #significance test for interaction
          mod_inter_ML <- rma(yi=yi, mods = ~time * group, vi=vi, data=es_group_dat, method = "ML") 
         mod_red_ML <- rma(yi=yi, mods = ~time + group, vi=vi, data=es_group_dat, method = "ML") 
  
          pval <- anova(mod_inter_ML, mod_red_ML)[6]
          pval_clean <- ifelse(pval$pval < 0.001, "< 0.001", round(pval$pval,4))
          
  # get OR of time from two models with each level as reference level
          inter_mod1 <- rma(yi=yi,mods = ~ time_div * group, vi=vi, data=es_group_dat, method = "REML")
          inter_mod2 <- rma(yi=yi,mods = ~ time_div * I(relevel(group, newref)) , vi=vi, data=es_group_dat, method = "REML")
          
          OR_data$OR <- c("", round(c(exp(inter_mod1$b[2]), exp(inter_mod2$b[2])),2))
          OR_data$lowCI <- c("", round(c(exp(inter_mod1$ci.lb[2]), exp(inter_mod2$ci.lb[2])),2))
          OR_data$upCI <- c("", round(c(exp(inter_mod1$ci.ub[2]), exp(inter_mod2$ci.ub[2])),2))
          OR_data$p <- c("", round(c(inter_mod1$pval[2], inter_mod2$pval[2]),4))
          OR_data$`p inter` <- c(pval_clean, "", "")
  
  return(OR_data)
}


########
### function to get OR of time and p-values per continent 
########

# NB = character string for the variable with sample size
# t = character string for the time variable
# y = character string for the variable with the number of events
# group = character string for the group variable
# newref1 = character string for one level of the group variable that is not the reference level
# newref2 = character string for another level of the group variable that is not the reference level
# var = character string for the name of the variable in the table

sub_OR_cont_fun <- function(NB, t, y, group, newref1, newref2, var, pred_mat, labX, labY, title, limY, ypos, vec){
  
  OR_data <- data.frame(Variable = c(var, rep("",3)), Continent = rep(NA,4), N = rep(NA,4), 
                                     OR = rep(NA,4), lowCI = rep(NA,4), upCI = rep(NA,4),
                                      "OR (95% CI)" = rep(NA,4),
                                     p = rep(NA, 4), "p inter" = rep(NA,4), check.names = FALSE)
  
  cont_dat <- na.omit(DATA_cont[,c(NB, y, t, group)])
  colnames(cont_dat)<-c("Ni","Ei","time","continent")
  cont_dat$time_div <- cont_dat$time/10
  cont_dat$prop <- (cont_dat$Ei/cont_dat$Ni) *100
  
  OR_data$Continent <- c("", c("Asia", "Europe", "North-America"))
  OR_data$N <- c("", table(cont_dat$continent))
  
  es_cont_dat <- escalc(
    xi= Ei,
    ni =Ni,
    data = cont_dat, 
    measure = "PLO",
    to="if0all")
  
  #significance test for interaction
  mod_inter_ML <- rma(yi=yi, mods = ~time * continent, vi=vi, data=es_cont_dat, method = "ML") 
  mod_red_ML <- rma(yi=yi, mods = ~time + continent, vi=vi, data=es_cont_dat, method = "ML") 
  
  pval <- anova(mod_inter_ML, mod_red_ML)[6]
  pval_clean <- ifelse(pval$pval < 0.001, "< 0.001", round(pval$pval,4))
  
  # get OR of time from two models with each level as reference level
  inter_mod1 <- rma(yi=yi,mods = ~ time_div * continent, vi=vi, data=es_cont_dat, method = "REML")
  inter_mod2 <- rma(yi=yi,mods = ~ time_div * I(relevel(continent, newref1)) , vi=vi, data=es_cont_dat, method = "REML")
  inter_mod3 <- rma(yi=yi,mods = ~ time_div * I(relevel(continent, newref2)) , vi=vi, data=es_cont_dat, method = "REML")
  
  
  OR_data$OR <- c("", round(c(exp(inter_mod1$b[2]), exp(inter_mod2$b[2]), exp(inter_mod3$b[2])),2))
  OR_data$lowCI <- c("", round(c(exp(inter_mod1$ci.lb[2]), exp(inter_mod2$ci.lb[2]), exp(inter_mod3$ci.lb[2])),2))
  OR_data$upCI <- c("", round(c(exp(inter_mod1$ci.ub[2]), exp(inter_mod2$ci.ub[2]), exp(inter_mod3$ci.ub[2])),2))
  OR_data$`OR (95% CI)` <- ifelse(OR_data$OR != "", paste0(OR_data$OR, " (", OR_data$lowCI," - ", OR_data$upCI, ")"), "")
  OR_data$p <- c("", round(c(inter_mod1$pval[2], inter_mod2$pval[2], inter_mod3$pval[2]),4))
  OR_data$`p inter` <- c(pval_clean, "", "", "")
  
  #Prediction
  inter_mod_full <- rma(yi=yi,mods = ~ time * continent, vi=vi, data=es_cont_dat, method = "REML")
  
  pred_inter <- as.data.frame(predict(inter_mod_full, newmods = pred_mat, transf=transf.ilogit, addx=T)) %>%
    mutate(continent = c(rep(levels(es_cont_dat$continent)[1], vec[1]), 
                     rep(levels(es_cont_dat$continent)[2], vec[2]),
                     rep(levels(es_cont_dat$continent)[3], vec[3])))
  
  #sample size specific data points
  es_cont_dat$size <-  2 + 20* ((1 / sqrt(es_cont_dat$vi)) - min((1 / sqrt(es_cont_dat$vi)),na.rm=T))/(max((1 / sqrt(es_cont_dat$vi)),na.rm=T) - min((1 / sqrt(es_cont_dat$vi)),na.rm=T))
  
  #Graph
  cont_plot <- ggplot(data=es_cont_dat, aes(x= time, y= prop, colour = continent)) +
    geom_line(data=pred_inter, aes(x= X.time, y= pred*100, colour = continent), size=2) +
    geom_line(data=pred_inter, aes(x= X.time, y= ci.lb*100, colour = continent), size=1,linetype = 2, alpha = 0.5) +
    geom_line(data=pred_inter, aes(x= X.time, y= ci.ub*100, colour = continent), size=1,linetype = 2, alpha = 0.5) +
    geom_jitter(size = es_cont_dat$size, shape=1, stroke=1, width = 0.3, alpha = 0.5) +
    scale_color_brewer(palette = "Dark2", name = "", labels = c("Asia", "Europe", "North-America")) +
    guides(size=FALSE) + ylab(labY) + xlab(labX) +
    theme_classic() +
    theme(axis.title.x = element_text( face="bold", size=30),
          axis.title.y = element_text( face="bold", size=30),
          axis.text.x = element_text(face="bold", size=30),
          axis.text.y = element_text(face="bold", size=30),
          plot.title = element_text(face = "bold",size = 30, hjust = 0.5),
          legend.text = element_text(size=20), legend.title = element_text(size = 22),
          legend.position = "top") + ggtitle(title) +
    ylim(limY)
  
  OR_cont1 <- paste0(OR_data$Continent[2]," OR (CI): ",OR_data$OR[2], " (", OR_data$lowCI[2],"; ", 
                   OR_data$upCI[2],")")
  OR_cont2 <- paste0(OR_data$Continent[3]," OR (CI): ",OR_data$OR[3], " (", OR_data$lowCI[3],"; ", 
                   OR_data$upCI[3],")")
  OR_cont3 <- paste0(OR_data$Continent[4]," OR (CI): ",OR_data$OR[4], " (", OR_data$lowCI[4],"; ", 
                     OR_data$upCI[4],")")
  
  OR_p_plot <- paste0(OR_cont1,"\n",OR_cont2, "\n", OR_cont3,"\n p interaction=", OR_data$`p inter`[1])
  
  grob <- grobTree(textGrob(OR_p_plot, x=0.02,  y=ypos, hjust=0,
                            gp=gpar(col="black", fontsize=16, fontface="bold")))
  
  cont_plot <- cont_plot + annotation_custom(grob)
  
  
  return(list(OR_data, cont_plot))
}


sub_nonlin_cont_fun <- function(NB, t, y, group, knots){
  
  cont_dat <- na.omit(DATA_cont[,c(NB, y, t, group)])
  colnames(cont_dat)<-c("Ni","Ei","time","continent")
  cont_dat$time_div <- cont_dat$time/10
  cont_dat$prop <- (cont_dat$Ei/cont_dat$Ni) *100
  
  es_cont_dat <- escalc(
    xi= Ei,
    ni =Ni,
    data = cont_dat, 
    measure = "PLO",
    to="if0all")
  
  #significance test for interaction
  mod_inter_ML <- rma(yi=yi, mods = ~ rcspline.eval(time, nk = knots,inclx = TRUE)*continent, 
                              vi=vi, data=es_cont_dat, method = "ML") 
  mod_red_ML <- rma(yi=yi, mods = ~ time * continent, 
                              vi=vi, data=es_cont_dat, method = "ML") 
  
  pval <- anova(mod_inter_ML, mod_red_ML)[6]
  pval_inter <- ifelse(pval$pval < 0.001, "< 0.001", round(pval$pval,4))
  
  return(pval_inter)
}