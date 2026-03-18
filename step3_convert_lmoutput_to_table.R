generate_table <- function(file){
  fh <- readLines(file)
  out <- cbind.data.frame("Parameter", "Estimate", "T value", "R²")
  
  factors = F
  for(i in 1:length(fh)){
    fh[i] <- gsub("\\s+", " ", fh[i])
    chr <- unlist(strsplit(fh[i], ""))
    col <- unlist(strsplit(fh[i], " "))
    if(length(chr) > 1){
      if(chr[1] == ":"){
        out <- rbind(out, c(col[2], "", "", ""))
        factors = F
      } else if(col[1] == "(Intercept)"){
        factors = T
      } else if(col[1] == "---" | col[1] == ""){
        factors = F
      } else if(factors){
        if(is.na(col[6])){
          col[6] <- ""
        }
        if(col[5] == "<" | col[5] == "<2e-16" | as.double(col[5]) <= 0.001){
          col[6] <- "***"
        } else if(as.double(col[5]) <= 0.01){
          col[6] <- "**"
        } else if(as.double(col[5]) <= 0.05){
          col[6] <- "*"
        }
        out <- rbind(out, c(col[1], paste0(round(as.double(col[2]), 3), "±", round(as.double(col[3]), 3)), paste0(round(as.double(col[4]), 3), " ", col[6]), ""))
      } else if(col[1] == "Multiple"){
        out[nrow(out), 4] <- round(as.double(gsub(",", "", col[3])), 3)
        if(out[nrow(out), 4] == 0){
          out[nrow(out), 4] <- 0.001
        }
      }
    }
  }
  return(out)
}

write.table(generate_table("lm_output.hexad.txt"), "lm_output.hexad.table.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
write.table(generate_table("lm_output.txt"), "Table2.S5.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
write.table(generate_table("lm_output-2.txt"), "Table2.S5.part2.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
write.table(generate_table("lm_output-3.txt"), "TableS3.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
write.table(generate_table("glm_output.txt"), "TableS6.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
write.table(generate_table("glm_output-2.txt"), "TableS6-part2.txt", quote = F, row.names = F, sep = "\t", fileEncoding="UTF-8")
