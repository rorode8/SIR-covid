# Curvas de comportamiento de COVID 19


[repositorio de github](https://github.com/rorode8/SIR-covid)




El comportamiento de epidemias puede ser simulado con modelos matemáticos o modelos estadísticos. En este estudio nos concentraremos en los modelos matemáticos. Específicamente usaremos un modelo **compartimental**. Los modelos de este tipo simplifican el problema al dividir la población total en compartimentos, asumiendo que cada individuo en un mismo compartimento tiene las mismas características.


  
# Identificación del Problema
## Objetivo


Desarrollar diferentes modelos matemáticos por medio de ODEs que nos permitan entender distintos aspectos de la pandemia, tales como medir la velocidad de propagación de los contagios, establecer el rol de algunos mecanismos en el proceso de propagación, evaluar medidas de mitigación y contención, hacer estimaciones de la capacidad hospitalaria, etc.


## Requisitos


Para simular el comportamiento de la epidemia usaremos el modelo **SIR **(**S**usceptible **I**nfectious **R**emoved)[1]. El modelo consiste de al menos tres compartimientos.



   -  **S**: Personas susceptibles. Estos se pueden pasar al grupo de **I**nfectados al tener contacto con uno de ellos. 
   -  **I**: Personas infectadas. Aquellos que son portadores de la enfermedad. 
   -  **R**: Personas con inmunidad o que murieron. Este grupo puede ser interpretado de varias formas, pero son aquellas que fueron infectadas y ya no pertenecen a ese grupo. 



Cabe mencionar que este modelo, a pesar de ser muy usado en epidemiología, es muy básico porque no considera nacimientos o muertes en el análisis, tampoco la edad de las personas o la situación o estructura social en la que se encuentran. 




Al final compararemos nuestras gráficas de variables continuas contra las variables discretas y veremos qué similitudes y diferencias se crean entre la realidad y este modelo matemático simple.




Para considerarse exitoso el modelo deberá asimilar los datos reales en al menos un escenario con datos previos contra los que se le compare.


# **Modelo Matemático**
  
### Modelo SIR


Usaremos el modelo **SIR** de acuerdo al siguiente sistema de ecuaciones:



<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\frac{\textrm{dS}}{\textrm{dt}}=-\beta&space;\cdot&space;\mathrm{S}\cdot&space;\mathrm{I}+\delta&space;\cdot&space;\mathrm{R}\\&space;\frac{\textrm{dI}}{\textrm{dt}}=\beta&space;\cdot&space;\mathrm{S}\cdot&space;\mathrm{I}-\gamma&space;\cdot&space;\mathrm{I}\\&space;\frac{\textrm{dR}}{\textrm{dt}}=\gamma&space;\cdot&space;\mathrm{I}-\delta&space;\cdot&space;\mathrm{R}&space;\end{array}"/>



Aquí <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta&space;\cdot&space;\mathrm{S}\cdot&space;\mathrm{I}"/> representa la cantidad de personas que son contagiadas, es decir, los que pasan de ser susceptibles a ser contagiados (infectados). Es por esto que la constante <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta"/> (beta) representa la tasa de infección. Se sabe que unas cepas de covid son más infecciosas que otras.




De la misma forma <img src="https://latex.codecogs.com/gif.latex?\inline&space;\delta"/> representa la pérdida de inmunidad. Por eso <img src="https://latex.codecogs.com/gif.latex?\inline&space;\delta&space;\cdot&space;\mathrm{R}"/> representa a aquellos que pasan del compartimiento de ***R****emoved*, a ser susceptibles otra vez, es decir que pueden volverse a infectar.




Por último <img src="https://latex.codecogs.com/gif.latex?\inline&space;\gamma"/>(gamma) representa la tasa de recuperación.




![image_0.png](README_images/image_0.png)




    **Figura 1: **Diagrama de flujo del modelo SIR


  
### Modelo SIRD


Podemos agregar más variantes al modelo base **SIR **para hacerlo más interesante y realista. En este caso agregamos un compartimento extra para los decesos. De esta forma, solo pasarán a **R** los recuperados y los muertos se pasarán al nuevo compartimento **D**




![image_1.png](README_images/image_1.png)




**    Figura 2: **Diagrama de flujo del modelo SIRD** **


  


Nuestro nuevo sistema de ecuaciones se ve de la siguiente manera.



<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\frac{\textrm{dS}}{\textrm{dt}}=-\beta&space;\cdot&space;\mathrm{S}\cdot&space;\mathrm{I}+\delta&space;\cdot&space;\mathrm{R}\\&space;\frac{\textrm{dI}}{\textrm{dt}}=\beta&space;\cdot&space;\mathrm{S}\cdot&space;\mathrm{I}-\gamma&space;\cdot&space;\mathrm{I}-\mu&space;\cdot&space;I\\&space;\frac{\textrm{dR}}{\textrm{dt}}=\gamma&space;\cdot&space;\mathrm{I}-\delta&space;\cdot&space;\mathrm{R}\\&space;\frac{\textrm{dD}}{\textrm{dt}}=\mu&space;\cdot&space;I&space;\end{array}"/>

## Creación del modelo con datos


Antes de implementar el modelo es importante mencionar las razones por las que el modelo no se verá como una curva real de covid con varios repuntes. Los coeficientes usados aquí son promedios obtenidos de los datos de [OWID](https://covid.ourworldindata.org/data/owid-covid-data.csv) para el Reino Unido (UK). Esto causa que no se tome en cuenta las diferentes medidas tomadas por el gobierno para controlar en cierto espacio de tiempo, solo refleja cómo se comporta en promedio.




Para la tasa de pérdida de inmunidad no existen los datos para calcularlo en para un país en específico, así que usaremos un valor promedio usado en [otros estudios](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7598795/).


  

```matlab:Code
% Model parameters for the UK
beta = 9.6905E-09; % rate of infection
gamma = 0.532271959; % rate of recovery 
delta = 1/60; % rate of immunity loss 100 days
mu = 0.085470461; % fatallity rate
N = 66650000; % Total UK population 
I0 = 10; % initial number of infected
T =  500; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)
fprintf('La tasa de reproducción R es de %.2f',N*beta/gamma)
```


```text:Output
La tasa de reproducción R es de 1.21
```


```matlab:Code
[t,S,I,R,D] = modelSIRD(beta,gamma,delta,mu,N,I0,T,dt);
```

# Pruebas y resultados

```matlab:Code
semilogy(t,[S,I,R,D])
hold on
%plot(t,I); 
%plot(t,R); 
%plot(t,D); 

%grid on;
xlabel('Days'); 
ylabel('Number of individuals');
legend('S','I','R','D','Location','southeast');
hold off
```


![figure_0.eps](README_images/figure_0.eps)


```matlab:Code

[v,i]=max(I);
d=i/4;
fprintf('El pico de casos activos es de %d y sucede %d días después de que empezó la pandemia',round(v),round(d))
```


```text:Output
El pico de casos activos es de 101819 y sucede 386 días después de que empezó la pandemia
```


```matlab:Code
[v,i]=max(D);
d=i/4;
fprintf('El total de muertes son %d al día %d desde que empezó la pandemia',round(v),round(d))
```


```text:Output
El total de muertes son 1412387 al día 500 desde que empezó la pandemia
```



La gráfica que resulta del modelo SIRD puede no distanciarse mucho de la gráfica real por las limitaciones que mencionamos antes. Es importante ver que el modelo solo nos muestra el rol que juega cada variable dentro de una determinada población. En el caso del Reino Unido, la tasa de mortalidad varía desde 2% hasta 14% dependiendo de la fecha.


## Comparación datos reales


Comparemos la acumulación de muertes de nuestra predicción contra los datos de OWID.



```matlab:Code
data= readtable('UK-cum-death.csv');
data([end-499:end],[2,4,5])
```

| |areaName|date|cumDailyNsoDeathsByDeathDate|
|:--:|:--:|:--:|:--:|
|1|'United Kingdom'|06/12/2021|153753|
|2|'United Kingdom'|06/11/2021|153742|
|3|'United Kingdom'|06/10/2021|153729|
|4|'United Kingdom'|06/09/2021|153714|
|5|'United Kingdom'|06/08/2021|153698|
|6|'United Kingdom'|06/07/2021|153683|
|7|'United Kingdom'|06/06/2021|153677|
|8|'United Kingdom'|06/05/2021|153664|
|9|'United Kingdom'|06/04/2021|153658|
|10|'United Kingdom'|06/03/2021|153646|
|11|'United Kingdom'|06/02/2021|153631|
|12|'United Kingdom'|06/01/2021|153619|
|13|'United Kingdom'|05/31/2021|153609|
|14|'United Kingdom'|05/30/2021|153598|
|15|'United Kingdom'|05/29/2021|153590|
|16|'United Kingdom'|05/28/2021|153578|
|17|'United Kingdom'|05/27/2021|153565|
|18|'United Kingdom'|05/26/2021|153557|
|19|'United Kingdom'|05/25/2021|153543|
|20|'United Kingdom'|05/24/2021|153534|
|21|'United Kingdom'|05/23/2021|153518|
|22|'United Kingdom'|05/22/2021|153503|
|23|'United Kingdom'|05/21/2021|153499|
|24|'United Kingdom'|05/20/2021|153483|
|25|'United Kingdom'|05/19/2021|153477|
|26|'United Kingdom'|05/18/2021|153473|
|27|'United Kingdom'|05/17/2021|153464|
|28|'United Kingdom'|05/16/2021|153452|
|29|'United Kingdom'|05/15/2021|153439|
|30|'United Kingdom'|05/14/2021|153425|
|31|'United Kingdom'|05/13/2021|153411|
|32|'United Kingdom'|05/12/2021|153397|
|33|'United Kingdom'|05/11/2021|153382|
|34|'United Kingdom'|05/10/2021|153362|
|35|'United Kingdom'|05/09/2021|153354|
|36|'United Kingdom'|05/08/2021|153344|
|37|'United Kingdom'|05/07/2021|153324|
|38|'United Kingdom'|05/06/2021|153307|
|39|'United Kingdom'|05/05/2021|153290|
|40|'United Kingdom'|05/04/2021|153269|
|41|'United Kingdom'|05/03/2021|153254|
|42|'United Kingdom'|05/02/2021|153238|
|43|'United Kingdom'|05/01/2021|153220|
|44|'United Kingdom'|04/30/2021|153205|
|45|'United Kingdom'|04/29/2021|153189|
|46|'United Kingdom'|04/28/2021|153172|
|47|'United Kingdom'|04/27/2021|153151|
|48|'United Kingdom'|04/26/2021|153126|
|49|'United Kingdom'|04/25/2021|153107|
|50|'United Kingdom'|04/24/2021|153082|
|51|'United Kingdom'|04/23/2021|153058|
|52|'United Kingdom'|04/22/2021|153029|
|53|'United Kingdom'|04/21/2021|152996|
|54|'United Kingdom'|04/20/2021|152961|
|55|'United Kingdom'|04/19/2021|152928|
|56|'United Kingdom'|04/18/2021|152905|
|57|'United Kingdom'|04/17/2021|152875|
|58|'United Kingdom'|04/16/2021|152835|
|59|'United Kingdom'|04/15/2021|152805|
|60|'United Kingdom'|04/14/2021|152764|
|61|'United Kingdom'|04/13/2021|152731|
|62|'United Kingdom'|04/12/2021|152685|
|63|'United Kingdom'|04/11/2021|152649|
|64|'United Kingdom'|04/10/2021|152610|
|65|'United Kingdom'|04/09/2021|152567|
|66|'United Kingdom'|04/08/2021|152529|
|67|'United Kingdom'|04/07/2021|152469|
|68|'United Kingdom'|04/06/2021|152425|
|69|'United Kingdom'|04/05/2021|152380|
|70|'United Kingdom'|04/04/2021|152338|
|71|'United Kingdom'|04/03/2021|152288|
|72|'United Kingdom'|04/02/2021|152248|
|73|'United Kingdom'|04/01/2021|152198|
|74|'United Kingdom'|03/31/2021|152137|
|75|'United Kingdom'|03/30/2021|152063|
|76|'United Kingdom'|03/29/2021|152009|
|77|'United Kingdom'|03/28/2021|151941|
|78|'United Kingdom'|03/27/2021|151894|
|79|'United Kingdom'|03/26/2021|151832|
|80|'United Kingdom'|03/25/2021|151760|
|81|'United Kingdom'|03/24/2021|151685|
|82|'United Kingdom'|03/23/2021|151616|
|83|'United Kingdom'|03/22/2021|151524|
|84|'United Kingdom'|03/21/2021|151422|
|85|'United Kingdom'|03/20/2021|151299|
|86|'United Kingdom'|03/19/2021|151200|
|87|'United Kingdom'|03/18/2021|151099|
|88|'United Kingdom'|03/17/2021|151000|
|89|'United Kingdom'|03/16/2021|150875|
|90|'United Kingdom'|03/15/2021|150763|
|91|'United Kingdom'|03/14/2021|150629|
|92|'United Kingdom'|03/13/2021|150497|
|93|'United Kingdom'|03/12/2021|150363|
|94|'United Kingdom'|03/11/2021|150229|
|95|'United Kingdom'|03/10/2021|150061|
|96|'United Kingdom'|03/09/2021|149902|
|97|'United Kingdom'|03/08/2021|149717|
|98|'United Kingdom'|03/07/2021|149524|
|99|'United Kingdom'|03/06/2021|149330|
|100|'United Kingdom'|03/05/2021|149123|


```matlab:Code
column1 = table2cell(data(end-499:end,5)); %500 dias
deaths = cell2mat(column1);
days = 500:-1:1;
semilogy(days,deaths,"r o")
hold on
semilogy(t,D,"k")
xline(346,"--");
title("prediction model vs real data",' ')
ylabel("number of people")
xlabel("days elapsed")
legend("predicted","real","346 days of provided data","Location","best")
hold off
```


![figure_1.eps](README_images/figure_1.eps)



Como podía esperarse, al tomar una muestra tan grande de datos como lo son 346 días y simplemente promediar las tasas de cambio incurrimos en un error bastante grande.




Este modelo **no** cumple todos los requisitos formulados previamente ya que:



   -  el modelo tiene demasiado error y no se asemeja a la curva de los datos reales 

## CASO ITALIA SIRD SIMPLE


El siguiente script usa datos del 1 de agosto de 2021 a 1 de noviembre de 2021 de Italia, recolectados igualmente del OWID.



```matlab:Code
% ITALY
% Model parameters for the UK
beta = 1.65869E-09; % rate of infection
%no mask
%beta = beta*1.1;

gamma = 0.076163777; % rate of recovery 
delta = 1/60; % rate of immunity loss 100 days
mu = 0.001166467; % fatallity rate
N = 66650000 % Total UK population N = S + I + R
```


```text:Output
N = 66650000
```


```matlab:Code
I0 = 6714; % initial number of infected = [diferencia de casos en 30 dias - (diferencia decesos en 30 dias)]
T =  180; % period of 300 days
dt = 1/12; % time interval of 6 hours (1/4 of a day)
D0=35146;
%R0=281707
R0=160686
```


```text:Output
R0 = 160686
```


```matlab:Code
fprintf('La tasa de reproducción R es de %.2f',N*beta/gamma)
```


```text:Output
La tasa de reproducción R es de 1.45
```


```matlab:Code

[t,S,I,R,D] = modelSIRD2(beta,gamma,delta,mu,N,I0,T,dt,R0,D0);
semilogy(t,[S,I,R,D])
legend('S','I','R','D','Location','best');
title("180 days italy prediction")
```


![figure_2.eps](README_images/figure_2.eps)


```matlab:Code
[v,i]=max(D);
d=i/12;
fprintf('El total de muertes son %d al día %d desde que empezó la pandemia \n',round(v),round(d))
```


```text:Output
El total de muertes son 116827 al día 180 desde que empezó la pandemia 
```


```matlab:Code
data = csvread('ITA-CUM-death.csv');
column1 = data(:,1);
days = 1:length(column1);
plot(t,D, "k","LineWidth",2)
hold on
plot(days,column1,"r o")
xline(90,"k--","data provided till day 90")
grid on
hold off
legend("predicted","real")
xlabel('Days since august 1st 2020'); 
ylabel("Death toll")
```


![figure_3.eps](README_images/figure_3.eps)



Podemos ver cómo el modelo simula mejor el comportamiento de las muertes acumuladas reales, sin embargo, las autoridades en Italia tomaron [medidas preventivas para reducir los contagios y las muertes](https://www.garda.com/crisis24/news-alerts/367781/italy-authorities-extend-covid-19-restrictions-until-september-7-update-38). 




Este modelo **sí **cumple todos los requisitos formulados previamente ya que:



   -  es un modelo matemático de ODEs 
   -  su comportamiento es parecido a los datos reales 
   -  los datos usados para calcular los coeficientes sólo son de 90 días 

  
## CASO ITALIA SIRD AVANZADO [2]


Implementaremos el modelo del siguiente diagrama [2]




![image_2.png](README_images/image_2.png)




    **Figura 3: **Diagrama SIRD avanzado




Este modelo debería tomar una forma más realista ya que toma en cuenta más variables, como puede apreciarse  en la figura 3, tenemos los siguientes compartimentos.



   -  Susceptibles 
   -  Expuestos 
   -  Asintomáticos 
   -  Infectados 
   -  confinado (guarda distancia en casa) 
   -  En cuarentena (infectado pero ya no infecta más gente) 
   -  Recuperados 
   -  Muertes 

  

```matlab:Code
%italy
N = 6.046*10^7;
alpha = 0.0194;
beta=7.567;
mu=2.278 * 10^-6;
nao=9.180  * 10^-7;
sigma=1.463^-3;
tao=1.109  * 10^-4;
epsilon=0.263;
gamma=0.021;
delta=	0.077;
lambda0 = 	0.157;
lambda1 = 0.025;
k0 = 0.779;
k1 = 0.061;
x0=0;
T = 180;
dt = 1/4; 
y0= [N; 1; 0; 1; 1; 0; 1; 0];
[t,S, A, C, E, I, R, Q, D] = modelSIRDAvanzado(alpha,beta,mu,nao,sigma,tao,epsilon,gamma,delta,lambda0,lambda1,k0,k1,N,y0,x0,dt,T);
plot(t,S,'-',t, A,'--',t, C,':',t, E,'-.',t, I,'-',t, R,'--',t, Q,':',t, D,'-')
%ylim([0 3*10^7])
ylabel("Población")
xlabel("dias")
title("Italia SIRD-avanzado desde el inicio de la epidemia")
legend('S', 'A', 'C', 'E', 'I', 'R', 'Q', 'D')
```


![figure_4.eps](README_images/figure_4.eps)


```matlab:Code
plot(t,S,'-',t, A,'--',t, C,':',t, E,'-.',t, I,'-',t, R,'--',t, Q,':',t, D,'y-')
ylim([0 1.5*10^7])
ylabel("Población")
xlabel("dias")
title("Italia SIRD-avanzado desde el inicio de la epidemia (zoom)",' ')
legend('S', 'A', 'C', 'E', 'I', 'R', 'Q', 'D')
```


![figure_5.eps](README_images/figure_5.eps)


```matlab:Code
%data = csvread('ITA-CUM-death-day0-180.csv');
T = readtable('ITA-CUM-death-day0-180.csv');
column1 = T{:,1};

%column1 = data(:,1);
days = 1:length(column1);
semilogy(t,D, "k","LineWidth",2)
hold on
semilogy(days,column1,"r o")
%grid on
hold off
title("Muertes en Italia desde el inicio de la pandemia")
legend("predicted","real",'Location','southeast')
xlabel('Days'); 
ylabel("Death toll");
```


![figure_6.eps](README_images/figure_6.eps)



Aunque los números de la predicción son mayores, el comportamiento es parecido, pues ambas gráficas se empiezan a aplanar a los 40-60 días.




Este modelo con estos datos en particular **no **cumple con los requisitos definidos ya que los números no son realistas. 


## Entendiendo la importancia de los coeficientes


Revisemos el caso del Reino Unido  otra vez con el modelo SIRD simple para observar incrementos en las variables de tasa de infección (beta) y tasa de recuperación (gamma).



```matlab:Code
% Model parameters for the UK
beta = 9.6905E-09; % rate of infection
gamma = 0.532271959; % rate of recovery 
delta = 1/60; % rate of immunity loss 100 days
mu = 0.085470461; % fatallity rate
N = 66650000; % Total UK population 
I0 = 10; % initial number of infected
T =  500; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)
[t,S,I,R,D] = modelSIRD(beta,gamma,delta,mu,N,I0,T,dt);
```

### BETA


Los datos tomados para este estudio son los promedios de esas variables desde que empezó el 1 de febrero de 2020 hasta el pico de casos el 9 de enero de 2021. Para estas últimas fechas el proceso de vacunación (que empezó en diciembre 13 de 2020) juega un papel importante en reducir el número de decesos y contagios.




Veamos ahora qué pasa si se relajan las medidas preventivas, es decir, un aumento en el coeficiente de infección **beta**



```matlab:Code

betaBase = beta
```


```text:Output
betaBase = 9.6905e-09
```


```matlab:Code
for i = 1.1:0.1:1.5
    figure;
    beta = betaBase * i;
    
    [t,S,I,R,D] = modelSIRD(beta,gamma,delta,mu,N,I0,T,dt);
    semilogy(t,[S,I,R,D])    
    
    hold on
    title("SIRD con beta al "+i*100+" %")
    xlabel('Days'); 
    ylabel('Number of individuals');
    legend('S','I','R','D','Location','southeast');
    hold off
    [v,i]=max(D);
    d=i/4;
    fprintf('el valor de beta es %d',beta)
    fprintf('La tasa de reproducción R es de %.2f',N*beta/gamma)
    fprintf('El total de muertes son %d al día %d desde que empezó la pandemia \n',round(v),round(d))
end
```


![figure_7.eps](README_images/figure_7.eps)


```text:Output
el valor de beta es 1.065955e-08
La tasa de reproducción R es de 1.33
El total de muertes son 5716279 al día 500 desde que empezó la pandemia 
```


![figure_8.eps](README_images/figure_8.eps)


```text:Output
el valor de beta es 1.162860e-08
La tasa de reproducción R es de 1.46
El total de muertes son 9269041 al día 500 desde que empezó la pandemia 
```


![figure_9.eps](README_images/figure_9.eps)


```text:Output
el valor de beta es 1.259765e-08
La tasa de reproducción R es de 1.58
El total de muertes son 12173548 al día 500 desde que empezó la pandemia 
```


![figure_10.eps](README_images/figure_10.eps)


```text:Output
el valor de beta es 1.356670e-08
La tasa de reproducción R es de 1.70
El total de muertes son 14777172 al día 500 desde que empezó la pandemia 
```


![figure_11.eps](README_images/figure_11.eps)


```text:Output
el valor de beta es 1.453575e-08
La tasa de reproducción R es de 1.82
El total de muertes son 16935175 al día 500 desde que empezó la pandemia 
```



Podemos ver cómo un incremento de apenas el 10% resulta en un gran incremento de decesos. En un escenario real, incrementa a la vez la tasa de mortalidad y decrementar por lo tanto la tasa de recuperación. Esto debido al sobrecupo de los hospitales.


### GAMMA

\hfill \break


```matlab:Code
beta = betaBase;
gammaBase = gamma;
for gamma = gammaBase:-0.11:gammaBase-0.11*4
    figure;
    
    [t,S,I,R,D] = modelSIRD(beta,gamma,delta,mu,N,I0,T,dt);
    semilogy(t,[S,I,R,D])    
    
    hold on
    title("SIRD con gamma = "+gamma)
    xlabel('Days'); 
    ylabel('Number of individuals');
    legend('S','I','R','D','Location','best');
    hold off
    [v,i]=max(D);
    d=i/4;
    fprintf('el valor de gamma es %d',gamma)
    fprintf('La tasa de reproducción R es de %.2f',N*beta/gamma)
    fprintf('El total de muertes son %d al día %d desde que empezó la pandemia \n',round(v),round(d))
end
```


![figure_12.eps](README_images/figure_12.eps)


```text:Output
el valor de gamma es 5.322720e-01
La tasa de reproducción R es de 1.21
El total de muertes son 1412387 al día 500 desde que empezó la pandemia 
```


![figure_13.eps](README_images/figure_13.eps)


```text:Output
el valor de gamma es 4.222720e-01
La tasa de reproducción R es de 1.53
El total de muertes son 10638058 al día 500 desde que empezó la pandemia 
```


![figure_14.eps](README_images/figure_14.eps)


```text:Output
el valor de gamma es 3.122720e-01
La tasa de reproducción R es de 2.07
El total de muertes son 21657465 al día 500 desde que empezó la pandemia 
```


![figure_15.eps](README_images/figure_15.eps)


```text:Output
el valor de gamma es 2.022720e-01
La tasa de reproducción R es de 3.19
El total de muertes son 34454982 al día 500 desde que empezó la pandemia 
```


![figure_16.eps](README_images/figure_16.eps)


```text:Output
el valor de gamma es 9.227196e-02
La tasa de reproducción R es de 7.00
El total de muertes son 48173566 al día 500 desde que empezó la pandemia 
```

# Conclusión


Entender este comportamiento de la pandemia causada por la enfermedad coronavirus nos permite ver las medidas sanitarias necesarias, estructurar medidas de mitigación y contención, hacer estimaciones de la capacidad hospitalaria, entre otros. Esto se debe a que en los modelos matemáticos usados podemos interactuar con las variables y así tratar de predecir el comportamiento de varios aspectos de la pandemia. En este caso es más fácil describir las variables por cómo cambian en lugar de por lo que representan. Es por esto que estos sistemas físicos son sistemas de ecuaciones diferenciales.




Como vimos en el caso de Inglaterra, al usar un promedio de un año para los datos, el modelo se distancia mucho de la realidad. A pesar de necesitar menos datos que los modelos estadísticos, los modelos de ecuaciones diferenciales ordinarias (como es el SIRD) dependen grandemente de los coeficientes.




En el primer caso de Italia con SIRD simple obtuvimos buenos resultados. Esto se debe a que los coeficientes tienen sentido para este marco de tiempo. Podríamos probablemente recortar los datos a tan solo 14 días y obtener una buena simulación de los siguientes 30 días si la situación no cambia. Es aquí donde vemos cómo estos modelos nos muestran un escenario realista o muy aproximado a la realidad. De esta forma, podemos reaccionar proactivamente antes de una catástrofe.




El último modelo está hecho para verse de una forma retrospectiva, difícilmente podremos calcular todos los coeficientes de forma previa. Pero podemos calcular diferentes escenarios jugando con estos coeficientes y ver la severidad de estos casos hipotéticos.




Finalmente, observamos cómo solamente el SIRD simple de Italia cumple con todos los requisitos. Al contrario del SIRD avanzado. Esto también nos muestra cómo los sistemas más complejos son muy delicados en cuanto a la exactitud de los parámetros (coeficientes) que se usan en el modelo. El modelo y su efectividad dependen más fuertemente de estos coeficientes que de su complejidad misma. 


# Referencias


[1] R. Beckley, C.Weatherspoon, M. Alexander , M. Chandler, A.Johnson , and G. S. Bhatt."Modeling epidemics with differential equations". Tennessee State University. [[link]](https://www.tnstate.edu/mathematics/mathreu/filesreu/GroupProjectSIR.pdf). visitado 5/7/2021




[2] [Ind. Eng. Chem. Res. 2021, 60, 11, 4251–4260](https://pubs.acs.org/doi/10.1021/acs.iecr.0c04754)


# Implementación

```matlab:Code
function [t,S,I,R] = modelSIR(beta,gamma,delta,N,I0,T,dt)
    [t,y] = runKut4(@(t,y) SIRfunc(t,y,beta,gamma,delta),0,dt,T,[N-I0;I0;0]);    
    S=y(:,1);
    I=y(:,2);
    R=y(:,3);
end

function [t,S,I,R,D] = modelSIRD(beta,gamma,delta,mu,N,I0,T,dt)
    [t,y] = runKut4(@(t,y) SIRDfunc(t,y,beta,gamma,delta,mu,N),0,dt,T,[N-I0;I0;0;0]);  
    S=y(:,1);
    I=y(:,2);
    R=y(:,3);
    D=y(:,4);
end

function [t,S,I,R,D] = modelSIRD2(beta,gamma,delta,mu,N,I0,T,dt,R0,D0)
    [t,y] = runKut4(@(t,y) SIRDfunc(t,y,beta,gamma,delta,mu,N),0,dt,T,[N-I0;I0;R0;D0]);  
    S=y(:,1);
    I=y(:,2);
    R=y(:,3);
    D=y(:,4);
end

function [t,S, A, C, E, I, R, Q, D] = modelSIRDAvanzado(alpha,beta,mu,nao,sigma,tao, ...
    epsilon,gamma,delta,lambda0,lambda1,k0,k1,N,y0,x0,dt,T)
    
    [t,y] = ode45(@(t,y) SIRDAfunct(t,y,alpha,beta,mu,nao,sigma,tao,epsilon ...
        ,gamma,delta,lambda0,lambda1,k0,k1,T,N),x0:dt:T,y0);
    
    S=y(:,1);
    A=y(:,2);
    C=y(:,3);
    E=y(:,4);
    I=y(:,5);
    R=y(:,6);
    Q=y(:,7);
    D=y(:,8);
    
end

function [xSol,ySol] = runKut4(F,x0,h,xStop,y0)
    
    reshape( y0, [],1 ); % column vector
    xSol = zeros(2,1); 
    ySol = zeros(length(y0),2);
    xSol(1) = x0; 
    ySol(:,1) = y0;
    i = 1;
    x=x0;
    y=y0;    
    while x < xStop
        i = i + 1;
        h = min(h,xStop - x);
        K1 = h*F(x,y);
        K2 = h*F(x + h/2,y + K1/2);
        K3 = h*F(x + h/2,y + K2/2);
        K4 = h*F(x+h,y + K3);
        y = y + (K1 + 2*K2 + 2*K3 + K4)/6;
        x = x + h;
        xSol(i) = x;
        ySol(:,i) = y; 
    end
    
    xSol = xSol';
    ySol=ySol';
end

function dydt = SIRfunc(t,y,beta,gamma,delta)
    dydt = zeros(3,1);
    %dydt(1)=dS
    %dydt(2)=dI
    %dydt(3)=dR
    dydt(1) = (-beta.*y(2).*y(1) + delta.*y(3));    
    dydt(2) = (beta.*y(2).*y(1) - gamma.*y(2));
    dydt(3) = (gamma.*y(2) - delta.*y(3));
end

function dydt = SIRDfunc(t,y,beta,gamma,delta,mu,N)
    dydt = zeros(4,1);
    %dydt(1)=dS
    %dydt(2)=dI
    %dydt(3)=dR
    dydt(1) = (-beta.*y(2).*y(1) + delta.*y(3));
    dydt(2) = (beta.*y(2).*y(1) - gamma.*y(2) -mu.*y(2));
    dydt(3) = (gamma.*y(2) - delta.*y(3));
    dydt(4) = mu.*y(2);
end

function dydt = SIRDAfunct(t,y,alpha,beta,mu,nao,sigma,tao,epsilon,gamma,delta,lambda0,lambda1,k0,k1,tf,N)
    dydt = zeros(8,1);
    %S=1, A=2, C=3, E=4, I=5, R=6, Q=7, D=8
    
    dydt(1) = -alpha*y(1) - beta*y(5)*y(6)./N - sigma*y(1)*y(2)./N -nao*y(1); %dS
    dydt(2) = - tao*y(2) + epsilon*y(4); %dA
    dydt(3) = (alpha*y(1)./N) -mu*y(3); %dC
    dydt(4) = -gamma*y(4) +beta*y(5)*y(1)./N +mu*y(3)+nao*y(1) +sigma*y(1)*y(2)./N -epsilon*y(4); %dE
    dydt(5) = tao*y(2)+gamma*y(4)-delta*y(5); %dI
    dydt(6) = (lambda0 + (t/tf)*(lambda1-lambda0))*y(7); %dR
    dydt(7) = delta*y(5)-((lambda0 + (t/tf)*(lambda1-lambda0))*y(7)) -((k0 + (t/tf)*(k1-k0))*y(7)); %dQ
    dydt(8) = (k0 + (t/tf)*(k1-k0))*y(7); %dD

end

```

