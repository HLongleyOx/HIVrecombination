##This script analyses linkage data derived from simulations to infer recombination according to the RATSLD method Romero and Feder 2024

#simulations recombination rate 
library(minpack.lm)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Set up input directory - modify this as needed
input_dir <- "linkage"  # Default input directory for linkage data
dir.create(input_dir, showWarnings = FALSE)

recom = c("5e-06", "7.5e-06", "1e-05","2.5e-05",  "5e-05")
df <- data.frame()

#function to fit data to curve as described in paper
f <- function(x, c0, c1, c2) {
  c0 + c1 * (1 - exp(-c2 * x))
}

recombination_estimate <-
  function(bootstrap,
           viralload_max,
           viralload_min,
           sel,
           timeframe_min,
           timeframe_max,
           subtype_choice,
           start,
           end,
           tech,
           max_distt){
    
    df_recombination <- data.frame()
    
    if (subtype_choice != "All") {
      bootstrap <- bootstrap %>% filter(subtype == subtype_choice)
    }
    
    bootstrap <- bootstrap %>% filter(
      siteA > start,
      siteB > start,
      siteA < end,
      siteB < end,
      viralload > viralload_min,
      viralload < viralload_max,
      last_days > timeframe_min,
      last_days < timeframe_max,
      freq_changeA<sel, 
      freq_changeB<sel
    ) 
    
    
    #bootstrap, pick one row per siteA, siteB, 
    for (i in 1:100){
        bootstrap_n = bootstrap
        boot_sample <- sample(1:(dim(bootstrap_n)[1]),(dim(bootstrap_n)[1]),replace=T)
        bootstrap_n <- bootstrap_n[boot_sample, ]
      
      linkage_recombination_nochange <-
        bootstrap_n  %>% filter(distt<max_distt) %>%
        dplyr::mutate(grp = round(distt /
                                    1) * 1) 
      x = linkage_recombination_nochange %>% group_by(grp) %>% summarise(mn = mean(ddratio), N = n())
      a = x$mn
      b = x$grp
      n = sqrt(x$N)
      
      tryCatch({
        fit <-
          nlsLM(
            a ~ f(b, c0, c1, c2),
            start = list(c0 = 0.2, c1 = 0.3, c2 = 0.0001),
            weights = n, control = list(maxiter = 500)
          )
        recombination = (coef(fit)[[2]] * coef(fit)[[3]])
        df = data.frame(
          rate = recombination,
          midpoint = start,
          N = dim(linkage_recombination_nochange)[1],
          run = i,
          viralload_min = viralload_min,
          viralload_max = viralload_max,
          timeframe_min = timeframe_min,
          timeframe_max = timeframe_max,
          sel = sel,
          tech = tech
        )
        df_recombination = rbind(df, df_recombination)
      }
      , error = function(err) {
        # Handle the error (optional) 
        print("error")
      })
    } 
    
    return(df_recombination)
  } 




for (file in recom){
  sel_linkage  <- read_csv(file.path(input_dir, paste0("linkage_", file, ".csv")), show_col_types = FALSE) %>%
    mutate(days=as.numeric(sub("_([0-9]+)_", "\\1", date)), siteA=as.numeric(siteA), siteB=as.numeric(siteB), pA=as.numeric(pA),  pB=as.numeric(pB), ddash=as.numeric(ddash),ld=as.numeric(ld)) %>% 
    drop_na() %>% group_by(ID, siteA, siteB, baseA, baseB) %>% filter(ld>=0.0001) 

  sel_first_timepoint <-  sel_linkage %>% filter(ddash>0.2) %>% group_by(ID,siteA, siteB, baseA, baseB) %>% filter(days==min(days)) %>% mutate(first_day = days) %>%
    ungroup() %>% select(siteA, siteB, first_day,ID) %>% right_join(sel_linkage) %>% filter(days>=first_day)  %>%
    group_by(ID,siteA, siteB, baseA, baseB) %>%
    arrange(days) %>%
    mutate(last_days = days-lag(days), dist=abs(siteA-siteB), ddash_lag=lag(ddash), lagPa=lag(pA), lagPb=lag(pB), distt=(last_days)*dist, ddratio = -1*log(ddash/ddash_lag), first_ld=first(ddash)) %>% 
    drop_na()
  
  

  for (d in seq(10000, 150000, by=10000)){
      sel <- sel_first_timepoint %>% drop_na() %>% filter(distt<d) %>%
        mutate(grp = round(distt /
                             1) * 1)
      
      x = sel %>% group_by(grp) %>% summarise(mn = mean(ddratio), N = n()) 
      a = x$mn
      b = x$grp
      n = sqrt(x$N)
      
      fit <-
        nlsLM(
          a ~ f(b, c0, c1, c2),
          start = list(c0 = 0.2, c1 = 0.3, c2 = 0.0001),
          weights = n
        )
      p=((coef(fit)[[2]] * coef(fit)[[3]]))
      print(p)
      df= rbind(df, data.frame(r=p, dist=d, true_recom=file))}

}


 df_50 <- df %>% filter(dist==50000, true_recom %in% c("1e-05", "2.5e-05","5e-05", "7.5e-06")) %>% mutate(x=r/as.numeric(true_recom))
cbPalette <- brewer.pal(4, "Set1")

sim_plot <- ggplot(df %>% filter(true_recom %in% c("1e-05", "2.5e-05","5e-05", "7.5e-06")) , aes(x=dist, y=r, group=as.factor(true_recom), colour=as.factor(true_recom))) + geom_line() + 
  theme_classic() + scale_y_log10(labels = scientific_10, breaks = c(1e-5, 3e-5, 5e-5,1e-4), limits=c(5*10^-6, 7.1*10^-5)) +
  ylab("Measured recombination rate\n(per site per generation)") + theme(axis.text = element_text(size=20), 
                                                         axis.title = element_text(size=20),
                                                         axis.title.x = element_text(margin = margin(t = 20)),
                                                         legend.text =  element_text(size=20),
                                                         legend.title =  element_text(size=20),
                                                         legend.position = "top") +
  xlab("Maximum time-scaled distance") +
  scale_colour_manual(values=cbPalette, name=expression(paste("True recombination rate (1e-5)")), labels=c("0.5", "0.75", "1","2.5","5")) +
  scale_x_continuous(breaks=c(50000, 100000), labels=c("50000","100000"))+
  geom_hline(data=df_50, aes(yintercept=as.numeric(true_recom), colour=true_recom), linetype=2, alpha=0.7) + guides(colour="none")


