ideal_rate_sampling_percentage = 5

resample_by_quantiles <- function(num_quantiles){
  regen_dt = as.data.table(do.call(rbind, transpose(approx(task_sums_1$tasks_percent, task_sums_1$hours_percent, seq(0, 100, 100 / num_quantiles)))))
  setnames(regen_dt,c("bins", "elapsed_times"))
  regen_dt$Pipeline = "SNV Regenotype"
  
  dup_dt = as.data.table(do.call(rbind, transpose(approx(task_sums_2$tasks_percent, task_sums_2$hours_percent, seq(0, 100, 100 / num_quantiles)))))
  setnames(dup_dt,c("bins", "elapsed_times"))
  dup_dt$Pipeline = "SV Dups"
  
  del_dt = as.data.table(do.call(rbind, transpose(approx(task_sums_3$tasks_percent, task_sums_3$hours_percent, seq(0, 100, 100 / num_quantiles)))))
  setnames(del_dt,c("bins", "elapsed_times"))
  del_dt$Pipeline = "SV Dels"
  
  disc_dt = as.data.table(do.call(rbind, transpose(approx(task_sums_4$tasks_percent, task_sums_4$hours_percent, seq(0, 100, 100 / num_quantiles)))))
  setnames(disc_dt,c("bins", "elapsed_times"))
  disc_dt$Pipeline = "SNV Discovery"
  
  return(rbind(regen_dt, dup_dt, del_dt, disc_dt))
}

butler_times = resample_by_quantiles(100 / ideal_rate_sampling_percentage)

butler_ideal_rates = butler_times[,.(bins, time_percent = elapsed_times - shift(elapsed_times, 1L, type="lag")), by=Pipeline][!is.na(time_percent)][,min(time_percent) / ideal_rate_sampling_percentage, by=Pipeline]
setkey(butler_ideal_rates, Pipeline)

butler_times = resample_by_quantiles(100)[,.(bins, time_percent = elapsed_times - shift(elapsed_times, 1L, type="lag")), by=Pipeline][!is.na(time_percent)]
setkey(butler_times, Pipeline)
butler_times = butler_times[butler_ideal_rates]
setnames(butler_times, "V1", "ideal_rate")
butler_times[, progress_rate := ideal_rate / time_percent, by=Pipeline]

ggplot(butler_times, aes(x=bins, y=progress_rate, color=Pipeline)) + 
  geom_smooth(span=0.3, se=FALSE) + 
  xlab("Percent Samples Completed") + 
  ylab("Actual Progress Rate / Target Progress Rate") +
  xlim(0,100) + 
  ylim(0,1.5) + 
  guides(col = guide_legend(title = "Pipeline")) + 
  theme(legend.background = element_rect(colour = "black"), legend.position=c(0.80,0.85), axis.title = element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title = element_text(size=18), axis.text = element_text(size=10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


main_ideal_rates = long_times[bins %% ideal_rate_sampling_percentage == 0, .(bins, time_percent = elapsed_times - shift(elapsed_times, 1L, type="lag")), by=Pipeline][!is.na(time_percent)][,min(time_percent) / ideal_rate_sampling_percentage, by=Pipeline]
setkey(main_ideal_rates, Pipeline)
setkey(long_times, Pipeline)
main_times = long_times[,.(bins, time_percent = elapsed_times - shift(elapsed_times, 1L, type="lag")), by=Pipeline][!is.na(time_percent)][main_ideal_rates]
setnames(main_times, "V1", "ideal_rate")
main_times[, progress_rate := ideal_rate / time_percent, by=Pipeline]

ggplot(main_times[Pipeline != 'DKFZ/EMBL'], aes(x=bins, y=progress_rate, color=Pipeline)) + 
  geom_smooth(span=0.3, se=FALSE) + 
  xlab("Percent Samples Completed") + 
  ylab("Actual Progress Rate / Target Progress Rate") +
  xlim(0,100) + 
  ylim(0,1.5) + 
  guides(col = guide_legend(title = "Pipeline")) + 
  theme(legend.background = element_rect(colour = "black"), legend.position=c(0.80,0.85), axis.title = element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title = element_text(size=18), axis.text = element_text(size=10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

mean_ratios = rbind(butler_times[,.(progress_rate = mean(progress_rate), pipeline="Butler"), by=bins],main_times[Pipeline != 'DKFZ/EMBL', .(progress_rate = mean(progress_rate), pipeline="PCAWG Core"), by=bins]) 
ggplot(mean_ratios, aes(x=bins, y=progress_rate, color=pipeline)) +
  geom_smooth(span=0.3, se=FALSE) +
  xlab("Percent Samples Completed") + 
  ylab("Mean Actual/Target Progress Rate") +
  xlim(0,100) + 
  ylim(0,1.5) + 
  guides(col = guide_legend(title = "Pipeline Manager")) + 
  theme(legend.background = element_rect(colour = "black"), legend.position=c(0.80,0.85), axis.title = element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title = element_text(size=18), axis.text = element_text(size=10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

  
  