library(data.table)
library(ggplot2)
library(scales)
library(fasttime)

elapsed_time_by_quantiles_2 <- function(data, num_quantiles, col_name){
  max_val = data[which(get(col_name) == max(get(col_name), na.rm=T)), col_name, with=F][1]
  bins = seq(0,num_quantiles)
  chunks =  bins * (max_val[[1]]/num_quantiles)
  chunk_indices = unlist(lapply(chunks, function(x) which.min(abs(unlist(data[,col_name, with=F]) - x))))
  max_time = as.numeric(data[which(get(col_name) == max(get(col_name), na.rm=T)), Date][1] - data[get(col_name) == 0, Date])
  elapsed_times = as.numeric(data[chunk_indices, Date] - data[get(col_name) == 0, Date]) * 100 / max_time
  return(data.table(bins, elapsed_times, col_name))
}

elapsed_time_by_quantiles <- function(data, num_quantiles, col_name){
  max_val = data[which(get(col_name) == max(get(col_name), na.rm=T)), col_name, with=F][1][[1]]
  max_val_index = which(data[,get(col_name)] == max_val)[1]
  min_val = 0
  max_time = as.numeric(data[which(get(col_name) == max(get(col_name), na.rm=T)), Date][1] - data[get(col_name) == 0, Date])
  val_indices = which(!is.na(data[1:max_val_index,get(col_name)]))
  percent_samples = data[val_indices, get(col_name)] * 100 / max_val
  percent_times = as.numeric(data[val_indices, Date] - data[get(col_name) == 0, Date]) * 100 / max_time
  
  intermediate = data.table(percent_samples, percent_times)

  result = as.data.table(do.call(rbind, transpose(approx(intermediate$percent_samples, intermediate$percent_times, seq(0, 100, 100 / num_quantiles)))))
  setnames(result,c("bins", "elapsed_times"))
  result[,Pipeline:=col_name]
  return(result)
}

#blah = elapsed_time_by_quantiles_2(prog, 100, "OxoG")

make_non_decreasing <- function(my_vec){
  val = 0
  
  for(i in seq_along(my_vec)){
    if(!(is.na(my_vec[i]) || is.na(val))){
      if(my_vec[i] < val){
        my_vec[i] = val
      }else{
        val = my_vec[i]
      }
    }
  }
  return(my_vec)
}

prog = fread("../data/pcawg_progress_sept_2016.csv")
prog$Date = as.Date(prog$Date, "%y-%m-%d")
prog$OxoG[586] = 0



col_list = c("BWA", "Sanger", "DKFZ/EMBL", "Broad", "OxoG")

for(el in col_list){
  prog[,el := make_non_decreasing(prog[,get(el),]), with=F]
}

prog[which(apply(prog, MARGIN =2, FUN=duplicated), arr.ind = T)] = NA

#col_list = c("OxoG")

times = lapply(col_list, function(x) elapsed_time_by_quantiles(prog, 100, x))
long_times = do.call(rbind, times)
setnames(long_times, "col_name", "Pipeline")
ggplot(long_times, aes(x=elapsed_times, y=bins, colour=Pipeline)) + 
  geom_line() + 
  xlab("Percent Runtime Elapsed") + 
  ylab("Percent Samples Complete") + 
  guides(col = guide_legend(title = "Pipeline")) + 
  theme(legend.background = element_rect(colour = "black"),
        legend.position=c(0.80,0.15), 
        axis.title = element_text(face="bold", size=16, family="Arial"), 
        legend.text=element_text(size=14, family="Arial"), 
        legend.title = element_text(size=18, family="Arial"), 
        axis.text = element_text(size=10)) + 
  geom_segment(x=0,y=0,xend=100,yend=100, linetype="dashed", colour="gray59") + 
  annotate("text", x=51, y=48, angle=45, label="Ideal Trajectory", colour="gray22",family="Arial") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



prog2 = melt(prog, id.vars=c("Date"), measure.vars = c("BWA", "Sanger", "DKFZ/EMBL", "Broad", "OxoG"), na.rm=T)


ggplot(prog2, aes(x=Date, y=value, colour=variable)) + 
  geom_line() + 
  scale_x_date("", labels=date_format("%b-%Y")) + ylab("# Donors")
