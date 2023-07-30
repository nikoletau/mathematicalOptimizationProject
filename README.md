# mathematicalOptimizationProject
Norwegian salmon supply chain

- model.py:contiene la costruzione del modello di ottimizzazione
- params.py raccoglie tutti i parametri dell'istanza che si sceglie di usare
- test.py: si crea l'istanza della classe Model a cui vengono assegnati i valori dei parametri; nel caso in cui il modello fosse infeasible, vengono restituiti i vincoli violati sullo standard output

- istanza 0: rappresenta il problema minimo, in cui vengono utilizzati solo 6 attori (uno per categoria); sono stati trovati i valori minimi (capacità, emissione) e massimi (distanza, domanda, consumo di carburante) per cui il problema minimo funziona 
- istanza 1: vengono utilizzati 8 attori (2 allevamenti di salmone, 1 macelleria, 1 impianto di elaborazione primaria, 1 impianto di elaborazione secondaria, 1 grossista e 2 rivenditori) e i periodi di tempo sono 3
- istanza 4: vengono utilizzati 20 attori (4 allevamenti di salmone, 3 macellerie, 2 impianti di elaborazione primaria, 2 impianti di elaborazione secondaria, 4 grossisti e 5 rivenditori) e i periodi di tempo sono 3
- istanza 7: vengono utilizzati 32 attori (10 allevamenti di salmone, 9 macellerie, 2 impianti di elaborazione primaria, 2 impianti di elaborazione secondaria, 4 grossisti e 5 rivenditori) e i periodi di tempo sono 14
- istanza 8: vengono utilizzati 46 attori (10 allevamenti di salmone, 9 macellerie, 8 impianti di elaborazione primaria, 10 impianti di elaborazione secondaria, 4 grossisti e 5 rivenditori) e i periodi di tempo sono 14
- istanza 9: vengono utilizzati 55 attori (10 allevamenti di salmone, 9 macellerie, 8 impianti di elaborazione primaria, 10 impianti di elaborazione secondaria, 8 grossisti e 10 rivenditori) e i periodi di tempo sono 14
- istanza 10: vengono utilizzati 95 attori (30 allevamenti di salmone, 9 macellerie, 8 impianti di elaborazione primaria, 10 impianti di elaborazione secondaria, 8 grossisti e 30 rivenditori) e i periodi di tempo sono 14
- istanza 11: vengono utilizzati 171 attori (60 allevamenti di salmone, 9 macellerie, 8 impianti di elaborazione primaria, 10 impianti di elaborazione secondaria, 24 grossisti e 60 rivenditori) e i periodi di tempo sono 14

Un periodo di tempo è un giorno. 

Per ogni istanza ci sono diversi .csv, in cui sono stati raggruppati i parametri a seconda degli indici. I .csv contengono i dati relativi ai collegamenti tra gli attori (capacità di trasporto, costo di trasporto, consumo di carburante (L/km) e distanza), mentre gli altri contengono la capacità di immagazzinamento ed eventuali costi di inventario ed elaborazione degli attori. 

Le istanze del test includono scenari in cui sono stati utilizzati diversi valori della domanda e sia le capacità del trasporto che quelle del magazzino sono state modificate quando il problema è risultato irrealizzabile o ha richiesto troppo tempo di computazione.

