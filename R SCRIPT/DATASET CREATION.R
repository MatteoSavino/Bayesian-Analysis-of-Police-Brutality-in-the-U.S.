##############################
#### DATASET CREATION  #######
##############################

{
  library(tidyverse)
  # for STAN
  library(rstan)
  library(rethinking)
  library(coda)
  # for regression
  library(car)
  library(rgl)
  # for excel files
  library(readxl)
  
  library(rgdal)
  
}

setwd('/Users/Francesco/Documents/POLICE_PROJECT')



#### WP DATASET & CENSUS ####

data <- read.csv('DATA/fatal_police_WP2.csv')      # fatal_police_kaggle/PoliceKillingsUS.csv   
usa <- read.csv('DATA/usa_population.csv')        # open USA_populaton.csv 


d <- data[which((data[,'armed'])!=''),]

Results<-c()
Results <- ifelse(d$armed !="unarmed" & d$race=="B", "BlackArmed",NA)
Results <- ifelse(d$armed =="unarmed" & d$race=="B", "BlackUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="H", "HispanicArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="H", "HispanicUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="W", "WhiteArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="W", "WhiteUnarmed",Results)
rm(d)
# or 
table(d$state,d$race) 

# We consider the states, not the counties!!!
table <- table(d$state,Results) 
dim(table)
table <- table[which(complete.cases(table[,'BlackArmed'])),] 
dim(table)

# Extract population size data for each state
# REM: Data from 2018:
usa <- read.csv('usa_population.csv',header=T) # open USA_populaton.csv 
dim(usa)

Nwhite <- usa$White
Nblack <- usa$Black
Nhispanic <- usa$Hispanic

# Extract shooting data for each state
UnarmedBlack <- as.numeric(table[,'BlackUnarmed'])
ArmedBlack <- as.numeric(table[,'BlackArmed'])
UnarmedHispanic <- as.numeric(table[,'HispanicUnarmed']	)
ArmedHispanic <- as.numeric(table[,'HispanicArmed'])
UnarmedWhite <- as.numeric(table[,'WhiteUnarmed']		)
ArmedWhite <- as.numeric(table[,'WhiteArmed'])



#### KAGGLE DATASET ####

# We have data by cities
geo_data1 <- read.csv('fatal_police_kaggle/MedianHouseholdIncome2015.csv')
geo_data2 <- read.csv('fatal_police_kaggle/PercentagePeopleBelowPovertyLevel.csv')
geo_data3 <- read.csv('fatal_police_kaggle/PercentOver25CompletedHighSchool.csv')
geo_data4 <- read.csv('fatal_police_kaggle/ShareRaceByCity.csv')

# Now we compute data by states:
Income <- geo_data1 %>%
  mutate(`Median Income`= as.numeric(`Median Income`)) %>%
  filter(is.na(`Median Income`)==FALSE) %>%
  group_by( `Geographic Area`) %>%
  summarize(income=sum(`Median Income`)/length(`Median Income`))

Poverty <- geo_data2 %>%
  mutate(poverty_rate= as.numeric(poverty_rate)) %>%
  filter(is.na(poverty_rate)==FALSE) %>%
  group_by( `Geographic Area`) %>%
  summarize(poverty=sum(poverty_rate)/length(poverty_rate))

Education <- geo_data3 %>%
  mutate(percent_completed_hs= as.numeric(percent_completed_hs)) %>%
  filter(is.na(percent_completed_hs)==FALSE) %>%
  group_by( `Geographic Area`) %>%
  summarize(education=sum(percent_completed_hs)/length(percent_completed_hs))

Ethnicity <- geo_data4 %>%
  mutate(share_black = as.numeric(share_black),
         share_white = as.numeric(share_white)) %>%
  filter(is.na(share_black)==FALSE) %>%
  filter(is.na(share_white)==FALSE) %>%
  group_by(`Geographic area`) %>%
  summarize(black=sum(share_black)/length(share_black),
            white=sum(share_white)/length(share_white))

Datageo <- cbind(Income,poverty=Poverty$poverty,education=Education$education,
                 black=Ethnicity$black,white=Ethnicity$white)

colnames(Datageo)[1] <- 'state'
head(Datageo)

rm(geo_data1,geo_data2,geo_data3,Income,Poverty,Education,Ethnicity)


########################## -------------------------------------------------------
#### ELECTION DATASET ####
##########################

election <- read_excel("DATA/usa_election2.xlsx") 
names(election)
#summary(election)
election <- election[,c(1:3,9,18,21,29:31,40:45)]
names(election)
colnames(election)[5] <- 'Population'
election$Population <- as.numeric(election$Population)

names(election)


### CALIFORNIA ------------------------------------------------

election <- election %>% filter(ST=='CA')
names(election)[4] <- 'Income'
head(election)
election <- election %>%
  mutate(Black = as.numeric(Black),
         Hispanic = as.numeric(Hispanic),
         White = as.numeric(White),
         Income = as.numeric(Income),
         Gini.Coefficient=as.numeric(Gini.Coefficient),
         Uninsured= as.numeric(Uninsured),
         Unemployment = as.numeric(Unemployment),
         Violent.crime = as.numeric(Violent.crime),
         Homicide.rate = as.numeric(Homicide.rate),
         Injury.deaths = as.numeric(Injury.deaths),
         Infant.mortality = as.numeric(Infant.mortality)) %>%
  mutate(Nblack = Population*Black,
         Nhispanic = Population*Hispanic,
         Nwhite = Population*White)

head(election)
save(election, file = 'california_election_data.RData')


# -------------------------------------------------------------
Gini <- election %>%
  filter(is.na(Gini.Coefficient)==0) %>%
  mutate(Gini.Coefficient = as.numeric(Gini.Coefficient)) %>%
  mutate(gini_tot = Gini.Coefficient/100*population) %>%
  group_by(ST) %>%
  summarise(Gini = sum(gini_tot)/sum(population)*100) 
dim(Gini)
Gini <- Gini[-52,2]


# We create a unique dataset:
dim(Gini)
dim(Datageo)
X <- cbind(Datageo,Gini)
save(X,file = 'DATA.RDAta')
rm(X)


########################## -------------------------------------------------------
####    CALIFORNIA    ####
##########################

# We compute data for California ------------------------------------------------------

# WashingtonPost Dataset and Population ----------------------------------------------

# <- read.csv('DATA/fatal_police_WP2.csv')      # fatal_police_kaggle/PoliceKillingsUS.csv   
usa <- read.csv('DATA/usa_population.csv')        # open USA_populaton.csv 

d <- data[which((data[,'armed'])!=''),]
d <- data %>% filter(state=='CA')
Results<-c()
Results <- ifelse(d$armed !="unarmed" & d$race=="B", "BlackArmed",NA)
Results <- ifelse(d$armed =="unarmed" & d$race=="B", "BlackUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="H", "HispanicArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="H", "HispanicUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="W", "WhiteArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="W", "WhiteUnarmed",Results)

# We consider the states, not the counties!!!
table <- table(d$city,Results) 
dim(table)
table <- table[which(complete.cases(table[,'BlackArmed'])),] 
dim(table)
head(table)

# KAGGLE DATASET ---------------------------------------------------

geo_data1 <- read.csv('fatal_police_kaggle/MedianHouseholdIncome2015.csv')
geo_data2 <- read.csv('fatal_police_kaggle/PercentagePeopleBelowPovertyLevel.csv')
geo_data3 <- read.csv('fatal_police_kaggle/PercentOver25CompletedHighSchool.csv')

geo_data1 <- geo_data1 %>% filter(Geographic.Area == 'CA')
geo_data2 <- geo_data2 %>% filter(Geographic.Area == 'CA')
geo_data3 <- geo_data3 %>% filter(Geographic.Area == 'CA')

geo_data <- cbind(geo_data1,poverty.rate=geo_data2[,3],HS.percentage=geo_data3[,3])
rm(geo_data1,geo_data2,geo_data3)

geo_data$City


# ELECTION DATASET ---------------------------------------------------------

election <- read_excel("DATA/usa_election.xlsx")
names(election)
#summary(election)
election <- election[,c(1,2,4,19,20,21,22,23,31,32,33,36,52,60:70)]
names(election)
colnames(election)[11] <- 'population'
election$population <- as.numeric(election$population)

election <- election %>% filter(ST=='CA') 

unique(election$city)


### LET US MERGE!!! ###

#ca_spdf <- readOGR('DATA/CA_Counties/CA_Counties_TIGER2016.shp')
CA_shp <- readOGR('/Users/Francesco/Documents/POLICE_PROJECT/DATA/CALIFORNIA/tl_2016_06_place.shp')
names(CA_shp@data)
head(CA_shp@data)

# We can create a link between our data!!!
unique(cali$city)
unique(CA_shp@data$NAME)
length(which(cali$city %in% CA_shp@data$NAME))

`%notin%` <- Negate(`%in%`)
missing <- setdiff(CA_shp@data$ST, geodata$ST)
good_row <- which(usashp@data$ST %notin% missing)
USA_spatial <- usashp[good_row, ]
USA_spatial@data <- data.frame(USA_spatial@data, geodata, victims, Nwhite=usa$White, Nblack=usa$Black, Nhispanic=usa$Hispanic,
                               Black.rate=usa$Black/usa$White*100, Hispanic.rate=usa$Hispanic/usa$White*100)  
USA_spatial@data$ST == geodata$ST
is(USA_spatial)





