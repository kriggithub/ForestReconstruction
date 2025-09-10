### Conclass function

data <- data.frame(
  col1 = c(2, 2, 3, 2, 5),
  col2 = c(4, 6, 5, 7, 2)
)



conclass <- function(status, decay){
  ifelse(status == 2 & decay %in%  1:2, 3, 
  ifelse(status == 2 & decay %in%  3:4, 4, 
  ifelse(status == 2 & decay %in%  5:6, 5, 
  ifelse(status == 3, 6, 
  ifelse(status == 4 & decay %in%  1:2, 3,
  ifelse(status == 4 & decay %in%  3:4, 4, 
  ifelse(status == 4 & decay %in%  5:6, 5, 
  ifelse(status == 4 & decay == 7, 8, 
  ifelse(status == 5, 7, NA)))))))))
}




data$col3 <- conclass(data$col1, data$col2)
