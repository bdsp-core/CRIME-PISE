crime.predict = function(
    t.PISE,       # PISE threshold 
    t.Death,      # Death threshold
    PISE.score,   # vector of predicted PISE score
    Death.score   # vector of predicted Death score
){
  df = data.frame(PISE.score, Death.score)
  colnames(df) = c("PISE_score", "Death_score")
  
  df.out = df %>%
    mutate(pred_label = case_when(
      PISE_score <= t.PISE & Death_score <= t.Death ~ "Event-Free",
      PISE_score > t.PISE & Death_score <= t.Death ~ "PISE",
      PISE_score <= t.PISE & Death_score > t.Death ~ "Death",
      PISE_score >= Death_score ~ "PISE",
      TRUE ~ "Death"
    ))
  
  # return the predicted label
  return(df.out$pred_label)
}


crime.evaluate = function(
    yr,   # at which year do you want to have the outcome evaluated? e.g., 1 year post-stroke
    t,    # time vector
    y,    # true label vector, with "PISE"/"Death"/"Censored"
    yhat, # predicted outcome, with "PISE"/"Death"/"Event-Free"
    PISE.score,   # predicted risk score for PISE
    Death.score   # predicted risk score for Death
){
  
  # recode actual outcome at yr years post-stroke
  yy = ifelse(t>=yr, "Event-Free", y)
  
  # when calculating classification result, only include patients with >yr year fu
  include = yy!="Censored"
  
  
  # time-dependent AUC (with competing risk) for PISE and Death
  ynum = recode(y, "Censored"=0, "PISE"=1, "Death"=2)
  troc.PISE = timeROC(
    t, 
    delta=ynum, 
    marker=PISE.score, 
    cause=1, 
    times=c(yr), 
    iid=T
  )
  troc.Death = timeROC(
    t, 
    delta=ynum, 
    marker=Death.score, 
    cause=2, 
    times=c(yr), 
    iid=T
  )
  
  # calculate classification evaluation metrics
  sen = sum(yy=="PISE" & yhat=="PISE" & include)/sum(yy=="PISE" & include)
  spe = sum(yy!="PISE" & yhat!="PISE" & include)/sum(yy!="PISE" & include)
  ppv = sum(yy=="PISE" & yhat=="PISE" & include)/sum(yhat=="PISE" & include)
  npv = sum(yy!="PISE" & yhat!="PISE" & include)/sum(yhat!="PISE" & include)
  n.pos = sum(yhat=="PISE" & include)
  n.pos.pise = sum(yy=="PISE" & yhat=="PISE" & include)
  n.pos.death = sum(yy=="Death" & yhat=="PISE" & include)
  n.pos.ef = sum(yy=="Event-Free" & yhat=="PISE" & include)
  
  # output eval metrics
  result = list(
    "AUC_PISE" = troc.PISE$AUC_2[[2]],
    "AUC_Death" = troc.Death$AUC_2[[2]],
    "sensitivity" = sen,
    "specificity" = spe,
    "ppv" = ppv,
    "npv" = npv,
    "n_pos" = n.pos,
    "n_pos_PISE" = n.pos.pise,
    "n_pos_Death" = n.pos.death,
    "n_pos_EventFree" = n.pos.ef
  )
  
  return(result)
}


make_cif_df = function(
    t,       # time vector
    y        # true label vector, with "PISE"/"Death"/"Censored"
){
  ynum = recode(y, "Censored"=0, "PISE"=1, "Death"=2)
  
  df = data.frame(t, ynum)
  colnames(df) = c("Time", "event")
  
  # fit a product limit function
  temp = prodlim(Hist(Time, event) ~ 1, data=df)
  
  # compute the probability of event at each time point
  df.pise = data.frame(c(0, temp$time), c(0, temp$cuminc$`1`), "PISE")
  df.death = data.frame(c(0, temp$time), c(0, temp$cuminc$`2`), "Death")
  colnames(df.pise) = c("Time", "Probability", "Group")
  colnames(df.death) = c("Time", "Probability", "Group")
  
  df.cif = rbind(df.pise, df.death)
  df.cif$Group = factor(df.cif$Group, level=c("PISE", "Death"))
  
  return(df.cif)
}


plot_individual_cif = function(
    pred.obj,   # predicted object
    lab_i,      # actual label event of this subject
    time_i,     # actual time of the event
    t,          # time vector of the reference population
    y           # true label vector of the reference population
    
){
  # get the training population risk reference
  df0 = make_cif_df(t, y)

  # get the model estimated CIF curve for every individual
  df1 = data.frame(c(0, pred.obj$time.interest), c(0, pred.obj$cif[,,1]), "PISE_i")
  df2 = data.frame(c(0, pred.obj$time.interest), c(0, pred.obj$cif[,,2]), "Death_i")
  colnames(df1) = c("Time", "Probability", "Group")
  colnames(df2) = c("Time", "Probability", "Group")
  
  # combine reference + individual data for visualization
  df.temp = rbind(df0, df1, df2)
  df.temp$Group = factor(df.temp$Group, levels=c("PISE_i", "PISE", "Death_i", "Death"))
  
  df.temp = subset(df.temp, Time<=5)
  g = ggplot(data=df.temp)+
    geom_step(aes(x=Time, y=Probability, color=Group))+
    ylim(c(0, 1))+
    theme_classic()+
    scale_color_manual(values=c("#e41a1c","#fbb4ae", "#377eb8", "#bdd7e7"))+
    scale_x_continuous(breaks = seq(0, 5, 1), limit=c(0, 5))
  
  if (lab_i=="PISE"){
    g = g+geom_vline(xintercept = time_i, color="#e41a1c", linetype="longdash")
  }else if (lab_i=="Death"){
    g = g+geom_vline(xintercept = time_i, color="#377eb8", linetype="longdash")
  }else{
    g = g+geom_vline(xintercept = time_i, color="grey", linetype="longdash")
  }
  
  return(g)
}


make_shap_df = function(
    mdl,   # model, i.e., SeQunatCR
    obs_i  # testing patient feature
){
  p_function1 <- function(model, data) predict(model, newdata = data)$predicted[,1]
  p_function2 <- function(model, data) predict(model, newdata = data)$predicted[,2]
  
  ive.result1 = individual_variable_effect(
    mdl, 
    data = mdl$xvar, 
    predict_function = p_function1, 
    new_observation = obs_i, 
    nsamples=50)
  
  ive.result2 = individual_variable_effect(
    mdl, 
    data = mdl$xvar, 
    predict_function = p_function2, 
    new_observation = obs_i,
    nsamples=50)
  
  return(list(
    "PISE" = ive.result1,
    "Death" = ive.result2
  ))
}


plot_individual_shap = function(
    shap.df   # the output of make_shap_df function
){
  
  # get the shap matrix of each outcome
  df.PISE = data.frame(
    shap.df$PISE$`_vname_`, 
    shap.df$PISE$`_attribution_`, 
    "PISE"
  )
  df.Death = data.frame(
    shap.df$Death$`_vname_`, 
    shap.df$Death$`_attribution_`, 
    "Death"
  )
  colnames(df.PISE) = c("Feature", "Contribution", "Outcome")
  colnames(df.Death) = c("Feature", "Contribution", "Outcome")
  
  # combine PISE and Death outcome for visualization
  df.temp = rbind(df.PISE, df.Death)
  df.temp$Outcome = factor(df.temp$Outcome, level=c("PISE", "Death"))
  df.temp$Feature = factor(df.temp$Feature, levels=shap.df$Death$`_vname_`)
  
  # make bar plot
  g = ggplot(df.temp, aes(x=Feature, y=Contribution, fill=Outcome))+
    geom_bar(stat="identity")+
    geom_hline(yintercept = 0)+
    coord_flip()+
    facet_wrap(.~Outcome)+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_fill_manual(values=c("#e41a1c","#377eb8"))
  
  return(g)
}
