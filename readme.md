# About this repo

This respository contains scripts for running simulations of pharmacogenomics data as described in the 2017 paper by Mistry and Chapman entitled "Simulation of cancer cell line pharmacogenomics data to optimise experimental design and analysis strategy" 

# Environment setup

Whilst this repo can be cloned and the scripts run on any computer, it is intended that the scripts be run on an Amazon Web Services instance for maximum reproducibility.  To set up such an environment follow the instructions below:

## Launch Amazon Instance with RStudio AMI

This tutorial assumes some familiarity with using AWS. Amazon's tutorials are very good, otherwise 
this [blog article](http://strimas.com/r/rstudio-cloud-1/) will help.  We used the [RStudio Server AMI's from Louis Aslett](http://www.louisaslett.com/RStudio_AMI/) from 25th May 2017 which comes with R 3.4.0 (specifically `ami-cbd3c6ad` on EU-West, Ireland).  

For now launch a `c4.2xlarge` instance with RStudio Server AMI.  Do this by clicking on the link from Louis Aslett's website, or by searching for the AMI identifier.  Under 'Configure Instance Details' the most imporant thing is to enable 'Auto-Assign Public IP' to allow you to access your instance, everything else should configure itself OK as default. Under 'Configure Security Group' ensure that ports 80 and 22 are open to all inbound and outbound traffic, this allows RStudio Server to work, and ssh access to your instance.  Now launch your instance, and either select or create a security key-pair (Amazon provide instructions on how to use these).  Monitor your instance as it is provisioned, and note the public IP address.  If it doesn't have a public IP address, you probably forgot to enable the 'Auto-Assign Public IP' option above.

## Connect to your AWS Instance

Once your instance has launched and configured, you should be able to connect to it in two ways:

- Using SSH with the following command: `ssh -i your-key-pair.pem ubuntu@your-instance-ip`.  You can download your keypair .pem file from AWS, and the IP of your instance will be visible in the instances dashboard.
- Using a web browser at `https://your-instance-ip/` and providing the username and password as `rstudio`.  This provides an RStudio Server session, which is like running RStudio on your desktop, except through a web browser where the computation is carried out on the cloud!

In principle from this point you could do everything you need in the rstudio environment, but it is helpful to configure the server a bit further at this point to install additional packages that are needed and make an AMI for future use.

## Configure your AWS Instance and make an AMI for future use

Connect to your server using ssh as described above and then switch to the root user using `sudo su -`.  Note that you have to be logged in as the `ubuntu` user to do this.  Now launch R and install the batchtools, tidyverse and devtools packages from CRAN using `install.packages`.  Finally install the `pgxsim` package using `devtools::install_github('chapmandu2/pgxsim')`.  This process can take 10-20 minutes.

Once everything is installed, reboot the AWS instance from the AWS Instances dashboard and connect again through the web browser.  Verify that everything is OK by running the following command:

```{r}
library(pgxsim)
example('sim_dose_response')
```

At this point you can create an image of your instance so that in future you can skip the above steps when you provision a new instance.  Do this by selecting Actions/Image/Create Image from the Instances dashboard.  You can also keep this particular instance for future use by Stopping the instance when you're not using it (this is advisible anyway!).  However, it is more expensive to store an instance than it is to store an AMI in the long term.

# Running the scripts

At this point you can make a copy of this git repo on your AWS instance by doing the following:

- Select File/New Project in RStudio Server
- Select Version Control/Git
- Paste the URL of this repo into the first box: `https://github.com/chapmandu2/pgx_simulation_scripts`
- Click Create

The simulation scripts need to be run interactively by stepping through them since they use the `batchtools` package to parallelise simulations. Initially it is advisible to change the `nreps` variable to a low number to ensure that everything works OK in a quick simulation. The scripts containing the word 'plots' can only be run once the main companion script is run, since it makes plots on the results of the simulations.  These scripts can also be run using  'File/Compile Report' which creates a pdf containing the plots, as well as the plots themselves.

For running the simulations in full we used the `c4.8xlarge` instance type which comes with 36 cores and costs \$1.87 per hour, for testing we used the smaller, cheaper instances such as `c4.2xlarge` with 8 cores at \$0.454 per hour.  On the 36 core instance with 200 simulation replicates, simulations 12 and 13 take about 30 minutes each, whereas simulation 11 takes around 5 hours.  Progress can be monitored either using Amazon's Cloudwatch or by connecting via ssh and running top.

# Running simulations on an HPC

The simulation screens need very little modification to run on an HPC or other environment since they use the `batchjobs` package.  The following two lines would have to be changed:

```
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = parallel::detectCores())
submitJobs(reg=reg, resources = list())
```

A different ClusterFunctions function would be specified, and the list object in resources would specify the number of CPU's etc to use.  Further information on configuring can be found on the batchjobs package website.


