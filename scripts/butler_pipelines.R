library(fasttime)
regen_task_instances = fread("../data/freebayes-discovery-regenotype-job-durations.tsv")
regen_task_instances[, end_date_days := cut(as.POSIXlt(regen_task_instances$end_date), "day")]
regen_cum_sums = regen_task_instances[order(end_date_days), .N, by="end_date_days"]
regen_cum_sums$cum_sum = cumsum(regen_cum_sums$N) / 24


discovery_task_instances = fread("../data/freebayes_discovery.tsv")

discovery_task_instances[, end_date_days := cut(as.POSIXlt(discovery_task_instances$end_date), "day")]
discovery_cum_sums = discovery_task_instances[order(end_date_days), .N, by="end_date_days"]
discovery_cum_sums$cum_sum = cumsum(discovery_cum_sums$N) / 24

del_task_instances = fread("../data/delly_deletion_regenotype.tsv")
del_task_instances[, end_date_days := cut(as.POSIXlt(del_task_instances$end_date), "day")]
del_cum_sums = del_task_instances[order(end_date_days), .N, by="end_date_days"]
del_cum_sums$cum_sum = cumsum(del_cum_sums$N)

dup_task_instances = fread("../data/delly_dup_regenotype.tsv")
dup_task_instances[, end_date_days := cut(as.POSIXlt(dup_task_instances$end_date), "day")]
dup_cum_sums = dup_task_instances[order(end_date_days), .N, by="end_date_days"]
dup_cum_sums$cum_sum = cumsum(dup_cum_sums$N)

get_task_sums <- function(task_instances){

  task_instances[,`:=`(execution_date = fastPOSIXct(execution_date),
                       start_date = fastPOSIXct(start_date),
                       end_date = fastPOSIXct(end_date))]
  min_start_date = min(task_instances$start_date)
  task_instances[, hours_elapsed := ceiling(as.double(difftime(end_date, min_start_date, units="hours")))]
  task_sums = task_instances[order(hours_elapsed), .N, by=.(hours_elapsed)][, task_sum := cumsum(N)]
  task_sums[, hours_percent := 100 * hours_elapsed / max(hours_elapsed)][, tasks_percent := 100 * task_sum / max(task_sum)]
  return(task_sums)
}

task_sums_1 = get_task_sums(regen_task_instances)

#Set initial point to 0.
task_sums_1 = rbind(task_sums_1, list(0,0,0,0,0))

task_sums_1$type = "SNV Regenotype"


task_sums_2 = get_task_sums(dup_task_instances)

#Set initial point to 0.
task_sums_2 = rbind(task_sums_2, list(0,0,0,0,0))

task_sums_2$type = "SV Duplications"


task_sums_3 = get_task_sums(del_task_instances)

#Set initial point to 0.
task_sums_3 = rbind(task_sums_3, list(0,0,0,0,0))

task_sums_3$type = "SV Deletions"

task_sums_4 = get_task_sums(discovery_task_instances)

#Set initial point to 0.
task_sums_4 = rbind(task_sums_4, list(0,0,0,0,0))

task_sums_4$type = "SNV Discovery"

task_sums = rbind(task_sums_1, task_sums_2, task_sums_3, task_sums_4)
ggplot(task_sums, aes(x=hours_percent, y=tasks_percent, color=type)) + 
  geom_line() + 
  xlab("Percent Runtime Elapsed") + 
  ylab("Percent Samples Complete") + 
  guides(col = guide_legend(title = "Pipeline")) + 
  theme(legend.background = element_rect(colour = "black"), legend.position=c(0.80,0.15), axis.title = element_text(face="bold", size=16), legend.text=element_text(size=14), legend.title = element_text(size=18), axis.text = element_text(size=10)) + 
  geom_segment(x=0,y=0,xend=100,yend=100, linetype="dashed", colour="gray59") + 
  annotate("text", x=45, y=40, angle=45, label="Ideal Trajectory", colour="gray22",family="Arial") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
