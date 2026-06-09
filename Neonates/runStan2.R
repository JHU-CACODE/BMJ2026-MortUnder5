
# Code to run a Stan as a background job
# This code uses the "st.input" list object of the Global environment

require(tidyverse)
require(rstan)
require(MCMCvis)


print(paste("Processing model:", st.input$name))

ptm0 = proc.time()

# First we prepare data suitable for the Stan model, using the information in the st.input list. Including the studies, deaths, lambda and other parameters

st.data <- list(
    # Studies
    nStudy = nrow(st.input$studies),
    K = ncol(st.input$studies) - 5,
    # matrix of vovariates
    Xmat = cbind(as.matrix(st.input$studies[,st.input$vxc]), scale(as.matrix(st.input$studies[,st.input$vxn]))),
    nNoshrink = st.input$nsv,
    first = st.input$studies$first,
    last = st.input$studies$last,
    # random effects
    nre = length(unique(st.input$studies$iso3)),
    reid = as.numeric(st.input$studies$iso3),
    Rnames=levels(st.input$studies$iso3), #names of Reffects
    # Deaths
    nCause = ncol(st.input$deaths) - 4,
    nRows = nrow(st.input$deaths),
    Missreport = as.matrix(st.input$deaths[,-c(1:4)]),
    nDeaths = st.input$deaths$n,
    # Parameters
    sd_betareg_noshrink = st.input$sdbeta,
    rsdlim = st.input$rsdlim,
    lambda = st.input$lambda,
    # Means and SD standardised covariates:
    xmeans = apply(st.input$studies[,st.input$vxn], 2, mean, na.rm=T),
    xsd = apply(st.input$studies[,st.input$vxn], 2, sd, na.rm=T)
  )

print(paste("Stan data created, now running model..."))

# Now fit the Stan model using the newly created st.data and other information from st.input
stanfit = rstan::sampling(st.input$model,
                          data = st.data,
                          pars = st.input$param,
                          chains = st.input$nchai, 
                          iter = st.input$niter, 
                          warmup = st.input$nwarm,
                          cores = st.input$cores,
                          control = list('adapt_delta' = .9),
                          seed = 1)

print(paste("Simulation finished, now computing summaries..."))

# if in st.input we stated that we only want the summary information from the model then substitute all the Stan output for only the summary table and add "_summary" to the name of the output object
if(st.input$summary==1){
  st.ouput <- MCMCsummary(stanfit,  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  st.input$name <- paste0(st.input$name,"_summary")
}else{
  st.output = stanfit
}

ptm1 = proc.time()
print(paste("Time taken (min):",(ptm1-ptm0)/60))

# Create object with outputs, inputs and Stan data
assign(st.input$name, list(st.output=st.output, st.data=st.data, st.input=st.input))
# Save the object in the specified folder
save(list=c(paste(st.input$name)), file=paste0(st.input$patho,"/",st.input$name,".RData"))






  
