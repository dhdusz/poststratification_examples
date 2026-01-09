/*designate library*/
libname kky "/knhanes";
libname bootlib "/knhanes_boot";


/* 1. Data preparation from KNHANES */ 

data revise ;
set kky.hn19_all;

/* (1) age category */
if age < 20 then age_cat = 0;
else if age >=20 and age <25 then age_cat = 2;
else if age >=25 and age <30 then age_cat = 2.5;
else if age >=30 and age <35 then age_cat = 3;
else if age >=35 and age <40 then age_cat = 3.5;
else if age >=40 and age <50 then age_cat = 4;
else if age >=50 and age <60 then age_cat = 5;
else if age >=60 and age <70 then age_cat = 6;
else if age >=70 and age <80 then age_cat = 7;
else if age =80  then age_cat = 8;

/* (2) obesity  */
if he_obe in (1,2,3) then obese = 0;
else if he_obe in (4,5,6) then obese = 1;

/* (3) past year alcohol consumption  */
if bd1 in (1,8) or bd1_11 = 1 then drink = 0;
else if bd1 = 2 and bd1_11 in (2,3,4,5,6) then drink = 1;

/* (4) lifetime cigarette : <5packs or >=5packs */
if bs1_1 in (1,3,8) then smoke = 0;
else if bs1_1 = 2 then smoke = 1;

/* (5) domain designation : note that age, sex, and region has no missing from the outset */

if age>=40 and drink ne . and smoke ne . then domain40 = 1;
else domain40 = 0;

if age<40  and age>=20 and drink ne . and smoke ne . then domain2030 = 1;
else domain2030 = 0;

if age>=40 and obese ne .  then domain_obese40 = 1;
else domain_obese40 = 0;

if age<40  and age>=20 and obese ne . then domain_obese2030 = 1;
else domain_obese2030 = 0;

keep ID sex age_cat obese drink smoke region kstrata psu wt_itvex domain40 domain2030 domain_obese40 domain_obese2030  ;

RUN;

/* 2. Get a point estimate of obesity for 2030 and 40+ */ 
PROC SURVEYFREQ DATA = revise nomcar ; strata kstrata; cluster psu; weight wt_itvex;
tables domain_obese40*obese /row ; 
tables domain_obese2030*obese / row;
ods output crosstabs= obesity ; 
RUN;

/* 3. Get marginal count of age_cat, sex, smoke, drink, region from the original KNHANES data*/
PROC SURVEYFREQ DATA = revise nomcar ; strata kstrata; cluster  psu; weight wt_itvex;
tables domain40*age_cat domain40*sex domain40*smoke domain40*drink domain40*region ; 
ods output crosstabs=test40 ; 
RUN;

PROC SURVEYFREQ DATA = revise nomcar ; strata kstrata; cluster  psu; weight wt_itvex;
tables domain2030*age_cat domain2030*sex domain2030*smoke domain2030*drink domain2030*region ; 
ods output crosstabs=test2030 ; 
RUN;

/* 3-1 Clean the marginal count results */
data temp2030_age ; 
	set test2030;
	if domain2030 = 1 and age_cat ne . and age_cat >=2 and age_cat < 4 ; 
	count = WgtFreq;
	label = "age";
	cat = age_cat;
	keep label cat count;
run;

data temp2030_sex ; 
	set test2030;
	if domain2030 = 1 and sex ne . ; 
	count = WgtFreq;
	label = 'sex' ; 
	cat = sex;
	keep label cat count;
run;

data temp2030_smoke ; 
	set test2030;
	if domain2030 = 1 and smoke ne . ; 
	count = WgtFreq;
	label = 'smk' ; 
	cat = smoke ; 
	keep  label cat count;
run;

data temp2030_drink ; 
	set test2030;
	if domain2030 = 1 and drink ne .  ; 
	count = WgtFreq;
	label = 'drk' ; 
	cat = drink ; 
	keep label cat count;
run;

data temp2030_region ; 
	set test2030;
	if domain2030 = 1 and region ne .  ; 
	count = WgtFreq;
	label = 'reg';
	cat = region; 
	keep label cat count;
run;

data temp2030_clean ; 
set temp2030_age temp2030_sex temp2030_smoke temp2030_drink temp2030_region; 
run;

data temp40_age ; 
	set test40;
	if domain40 = 1 and age_cat ne . and age_cat >= 4 ; 
	count = WgtFreq;
	label = "age";
	cat = age_cat;
	keep label cat count;
run;

data temp40_sex ; 
	set test40;
	if domain40 = 1 and sex ne . ; 
	count = WgtFreq;
	label = "sex";
	cat = sex;
	keep label cat count;
run;

data temp40_smoke ; 
	set test40;
	if domain40 = 1 and smoke ne . ; 
	count = WgtFreq;
	label = "smk";
	cat = smoke;
	keep label cat  count;
run;

data temp40_drink ; 
	set test40;
	if domain40 = 1 and drink ne .  ; 
	count = WgtFreq;
	label = "drk";
	cat = drink;
	keep label cat count;
run;

data temp40_region ; 
	set test40;
	if domain40 = 1 and region ne .  ; 
	count = WgtFreq;
	label = "reg";
	cat = region;
	keep label cat count;
run;

data temp40_clean ; 
set temp40_age temp40_sex temp40_smoke temp40_drink temp40_region; 
run;


proc export data=temp2030_clean
outfile = "/knhanes/knhanes_marginal_2030.csv"
dbms = csv replace ;
run;

proc export data=temp40_clean
outfile = "/knhanes/knhanes_marginal_40.csv"
dbms = csv replace ;
run;





/* 4. cluster bootstrap */

%macro rw_bootstrap(data=knhanes, nboot=1000, seed=12345, out_prefix=boot_);
option nonotes;
%do b = 1 %to &nboot;

   /* Step 1: Get list of unique PSUs by stratum */
   proc sql;
      create table psu_list as
      select distinct kstrata, psu
      from &data;
   quit;

   /* Step 2: Count PSUs per stratum */
   proc sql;
      create table strata_info as
      select kstrata,
             count(*) as n_psu,
             calculated n_psu / (calculated n_psu - 1) as rescale_wt
      from psu_list
      group by kstrata;
   quit;

   /* Step 3: Sample n-1 PSUs per stratum */
   data _null_;
      set strata_info;
      n_sample = n_psu - 1;
      call symputx('n_'||strip(put(_n_,best.)), n_sample);
      call symputx('strata_'||strip(put(_n_,best.)), kstrata);
      call symputx('nstrata', _n_);
   run;

   /* Step 4: Sample within each stratum separately */
proc sql;
      create table sampled_psus as
      select * from psu_list
      where 0=1;  /* Gets structure (column types) but zero rows */
   quit;

   %do s = 1 %to &nstrata;
      
      proc surveyselect data=psu_list(where=(kstrata=&&strata_&s))
                        method=urs
                        n=&&n_&s
                        seed=%eval(&seed + &b * 1000 + &s)
                        out=temp_psu(drop = NumberHits)
                        outhits
                        noprint;
      run;

   proc append base=sampled_psus data=temp_psu force;
   run;

   %end;

   /* Step 5: Merge selected PSUs with rescaling weights */
   proc sql;
      create table &out_prefix.&b as
      select &b as boot_iter, a.*,
             c.rescale_wt,
             a.wt_itvex * c.rescale_wt as bootstrap_weight
      from &data as a
      inner join sampled_psus as b
         on a.kstrata = b.kstrata and a.psu = b.psu
      inner join strata_info as c
         on a.kstrata = c.kstrata;
   quit;

   /* Clean up */
   proc datasets library=work nolist;
      delete psu_list strata_info sampled_psus temp_psu;
   quit;

%end;
option notes;
%mend;

%rw_bootstrap(data=revise, nboot = 1000, seed = 1);

/*4-1 save bootstrap results to bootlib library */
proc datasets library=work nolist;
   copy out=bootlib;
   select boot_1-boot_1000;  /* Range selection */
quit;



/* 5. Get marginal count of age_cat, sex, smoke, drink, region from bootstrapped KNHANES data  */

data counts_2030;
input boot_iteration label $3. cat count;
datalines ; ;
run; 

data counts_40;
input boot_iteration label $3. cat count;
datalines ; ;
run; 


%macro boot_counts(start = 1, stop=1000, out_prefix = boot_);
option nonotes ;

%do b = &start %to &stop;
ods select none;

data bootlib.&out_prefix.&b ; 
set bootlib.&out_prefix.&b ;
if age_cat >=4 and drink ne . and smoke ne . and obese ne . then domain_acc = 1;
else domain_acc = 0;
run;


proc surveyfreq data=bootlib.&out_prefix.&b nomcar ;
strata kstrata ; 
cluster psu ; 
weight bootstrap_weight ; 
tables domain2030*age_cat domain2030*sex domain2030*smoke domain2030*drink domain2030*region;  
ods output crosstabs = temp2030 ; 
run;

proc surveyfreq data=bootlib.&out_prefix.&b nomcar ;
strata kstrata ; 
cluster psu ; 
weight bootstrap_weight ; 
tables domain40*age_cat domain40*sex domain40*smoke domain40*drink domain40*region;  
ods output crosstabs = temp40 ; 
run;


ods select all;

data temp2030_age ; 
	set temp2030;
	if domain2030 = 1 and age_cat ne . and age_cat >=2 and age_cat < 4 ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "age";
	cat = age_cat;
	keep boot_iteration label cat count;
run;

data temp2030_sex ; 
	set temp2030;
	if domain2030 = 1 and sex ne . ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = 'sex' ; 
	cat = sex;
	keep boot_iteration label cat count;
run;

data temp2030_smoke ; 
	set temp2030;
	if domain2030 = 1 and smoke ne . ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = 'smk' ; 
	cat = smoke ; 
	keep boot_iteration label cat count;
run;

data temp2030_drink ; 
	set temp2030;
	if domain2030 = 1 and drink ne .  ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = 'drk' ; 
	cat = drink ; 
	keep boot_iteration label cat count;
run;

data temp2030_region ; 
	set temp2030;
	if domain2030 = 1 and region ne .  ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = 'reg';
	cat = region; 
	keep boot_iteration label cat count;
run;

data temp2030_clean_&b ; 
set temp2030_age temp2030_sex temp2030_smoke temp2030_drink temp2030_region; 
run;

data temp40_age ; 
	set temp40;
	if domain40 = 1 and age_cat ne . and age_cat >= 4 ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "age";
	cat = age_cat;
	keep boot_iteration label cat count;
run;

data temp40_sex ; 
	set temp40;
	if domain40 = 1 and sex ne . ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "sex";
	cat = sex;
	keep boot_iteration label cat count;
run;

data temp40_smoke ; 
	set temp40;
	if domain40 = 1 and smoke ne . ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "smk";
	cat = smoke;
	keep boot_iteration label cat  count;
run;

data temp40_drink ; 
	set temp40;
	if domain40 = 1 and drink ne .  ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "drk";
	cat = drink;
	keep boot_iteration label cat count;
run;

data temp40_region ; 
	set temp40;
	if domain40 = 1 and region ne .  ; 
	boot_iteration = &b;
	count = WgtFreq;
	label = "reg";
	cat = region;
	keep boot_iteration label cat count;
run;

data temp40_clean_&b ; 
set temp40_age temp40_sex temp40_smoke temp40_drink temp40_region; 
run;

proc append base = counts_2030 data = temp2030_clean force ; 
run; 

proc append base = counts_40 data = temp40_clean force ; 
run; 


%end;
option notes ;
%mend;

%boot_counts(start = 1, stop=1000);


proc export data=counts_2030
outfile = "/knhanes/knhanes_boot_marginal_2030.csv"
dbms = csv replace ;
run;

proc export data=counts_40
outfile = "/knhanes/knhanes_boot_marginal_40.csv"
dbms = csv replace ;
run;

