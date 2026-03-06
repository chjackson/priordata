
#' @export
fudge01 <- function(x, eps=1e-08){pmin(pmax(x, eps), 1 - eps)} # ?? used? 

# or could use S3 methods
pdata_is_normal <- function(pdata){
  is.list(pdata) && all(c("ybar","se") %in% names(pdata))
}
pdata_is_binomial <- function(pdata){
  is.list(pdata) && all(c("y","n") %in% names(pdata))
}




## todo consistent type for pdata is it list or df?
## generous in accept, strict in return
## df when we can conceive multiple rows

## todo work out how to structure when we have access to our own model
## sam <- model_fn(pdata)




