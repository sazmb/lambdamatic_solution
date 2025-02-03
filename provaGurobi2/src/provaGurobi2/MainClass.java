package provaGurobi2;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import gurobi.*;



	import gurobi.GRB;
	import gurobi.GRB.DoubleAttr;
	import gurobi.GRB.StringAttr;
	import gurobi.GRBEnv;
	import gurobi.GRBException;
	import gurobi.GRBLinExpr;
	import gurobi.GRBModel;
	import gurobi.GRBVar;

	public class MainClass 
	{
		static int n;
		static int m;
		static int k;
		static double c;
		static int h;
		static int x;
		static int y;
		private static boolean leggendoMatricedij;
		static int[] produzione;
		static int[] domanda_clienti;
		static int[][] distanze;
		static int [] alfa_j;
		

		public static void main(String[] args) {
			
			fileToVariable( );//caricare da file.txt
			try 
			{
			//creazione delle variabili e riempimento dei modelli	
				GRBEnv env = new GRBEnv("trasporto.log");
				impostaParametri(env);
				GRBModel model = new GRBModel(env);
				GRBModel modelIntero = new GRBModel(env);
				GRBModel rilassato = new GRBModel(env);
                
				GRBVar[][] xij = aggiungiVariabili(model, produzione, domanda_clienti);
				aggiungiFunzioneObiettivo(model, xij, distanze, c);
				GRBConstr[] constrProd=aggiungiVincoliProduzione(model, xij, produzione);
				GRBConstr[][] constrRag=aggiungiVincoliRaggio(model, xij, distanze,k);
				GRBConstr[] constrDom=aggiungiVincoliDomanda(model, xij, domanda_clienti);
				
				
				
				GRBVar[][] xijInt = aggiungiVariabili(modelIntero, produzione, domanda_clienti);
				GRBVar[] yi = aggiungiVariabiliIntere(modelIntero, alfa_j);
				aggiungiFunzioneObiettivoIntera(modelIntero, yi, alfa_j);
				aggiungiVincoliProduzioneIntera(modelIntero,xijInt, yi, produzione);
				aggiungiVincoliRaggio(modelIntero, xijInt, distanze,k);
				aggiungiVincoliDomanda(modelIntero, xijInt, domanda_clienti);
		        
				System.out.println("GRUPP0<13>");
				System.out.println("<Zambelli>");
				rispondiPrimaParte(model);
				rispondiSecondaParte(model,  distanze, h, xij, constrProd, constrRag, constrDom );
				rispondiTerzaParte(modelIntero, yi, xijInt, rilassato);
				
			}
			catch (GRBException e) 
			{
				e.printStackTrace();
			}	
			
		}
		
// start region PRIMA PARTE=======================================================================================================================================================
		private static void rispondiPrimaParte(GRBModel model) throws GRBException{
			System.out.println("\nQUESITO I");
			stampaSol(model);
			isDegenere(model);
			isMultipla(model);
		}
		
	   /*
		*  Verifica se la soluzione è degenere
	    *   conta il numero di variabili di base con valore uguale a zero
	    *   se anche solo una è maggiore stampa che è degenere
	    *   effettua anche un controllo anche nel caso limite in cui la soluzione ottima
	    *   si abbia con tutte le di slack in base 
	    *   
	    *  
	    */
		private static void isDegenere(GRBModel model)throws GRBException {
			 int count = 0;// numero di variabili in base con val==0
			    for (GRBVar var : model.getVars()) 
			    	 if (var.get(GRB.IntAttr.VBasis)==0 && roundToFourthDecimal(var.get(DoubleAttr.X))== 0) 
			    	        	 count++;
			    if (count != 0) 
			      System.out.println("Degenere:si");
			     else {
						for(GRBConstr constr: model.getConstrs())
							if(constr.get(GRB.IntAttr.VBasis)==GRB.BASIC && roundToFourthDecimal(constr.get(DoubleAttr.Slack))==0)
								count++;
						if(count!=0)
							System.out.println("Degenere:si");
						 else 
			    	 System.out.println("Degenere:no");
						 }
		}
		
		/*
		 * risolvo e stampo il valore della funzione obb. del problema
		 * ciclo che itera su tutte le var del model e ne stampa il valore
		 */	
		private static void stampaSol(GRBModel model) throws GRBException{
	model.set(GRB.IntParam.OutputFlag, 0);
	model.optimize();
	System.out.println("funzione obbiettivo=" +roundToFourthDecimal(model.get(GRB.DoubleAttr.ObjVal))+">");
	for(GRBVar var : model.getVars())
		System.out.println("<"+var.get(StringAttr.VarName)+">=<"+ roundToFourthDecimal(var.get(DoubleAttr.X))+">");
	
		
		}
	/*
	 * questo metodo va prima a iterare su tutte le variabili del modello se una variabile è
	 * in base con CR==0 aumenta il contatore 
	 * se il contatore !=0 stampa la risposta
	 * se no continua a indagare sui vincoli(per variabili di slack) e aumenta il contatore se questa è in base
	 *con CR( dato dall'attributo Pi del vincolo)
	 *se il contatore !=0 è multipla
	 *se no no
	 */
		private static void isMultipla(GRBModel model) throws GRBException{
			int count=0;
			for(GRBVar var:model.getVars())
				if(var.get(GRB.IntAttr.VBasis)!=GRB.BASIC && roundToFourthDecimal(var.get(GRB.DoubleAttr.RC))==0)
					count++;
			if(count!=0)
				System.out.println("Multipla:si");
			else {
				for(GRBConstr constr: model.getConstrs())
					if(constr.get(GRB.IntAttr.VBasis)!=GRB.BASIC && roundToFourthDecimal(constr.get(DoubleAttr.Pi))==0)
						count++;
				if(count!=0)
					System.out.println("Multipla:si");
				else
						System.out.println("Multipla:no");
				
			}
		}
//end region

//start region SECONDA PARTE==========================================================================================================
		private static void rispondiSecondaParte(GRBModel model, int [][] distanze, int h, GRBVar [][] xij, GRBConstr[] constrProd, GRBConstr[][] constrRag, GRBConstr[] constrDom )throws GRBException{
			System.out.println("\nQUESITO 2:");
			trovaVincoliNonAttivi(model, constrProd,  constrRag, constrDom);
			calcolaVariazioneK(model, distanze,xij);
			calcolaVariabiliDuale(model);
			calcolaVariazioneSh(model, h);
			calcolaVariazioneDxy(model, x ,y);
		}
		/*
		 * partendo da zero tenta di risolvere il problema per ogni istanza di p fino a che il model 
		 * non diventa risolvibile, ad ogni passo incrementa p di uno
		 * 
		 */
		private static void calcolaVariazioneK(GRBModel model, int [][] distanze, GRBVar [][] xij)throws GRBException{
     
			
	          int p=0;
	          do {
	        	  modificaVincoliRaggio(model, xij, distanze, p);//
	        	  model.update();
	        	  model.set(GRB.IntParam.OutputFlag, 0);
	        	  model.optimize();
	        	  p++;
	        	  
	          }while ((model.get(GRB.IntAttr.Status) == GRB.Status.INF_OR_UNBD ||
	        		    model.get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE ||
	        		    model.get(GRB.IntAttr.Status) == GRB.Status.UNBOUNDED));
	          p--;
	          System.out.println("intervallo k = [-INF,"+p+") ");

			
		}
		
		/*
		 * prima itera su ogni vincolo di <= del problema e verifica che il valore delle slack sia >0
		 * stessa cosa per i vincoli di >=ma stavolta controllando che sia <0
		 */
		 private static void trovaVincoliNonAttivi(GRBModel model,GRBConstr[] constrProd, GRBConstr[][] constrRag, GRBConstr[] constrDom)throws GRBException{
			System.out.printf("lista vincoli non attivi=");
			 
		
			for (GRBConstr c : constrProd) {
			    if (c.get(GRB.DoubleAttr.Slack) > 0.0) {
			    	  System.out.printf("[<" + c.get(GRB.StringAttr.ConstrName) + ">]");
			    }
			}
			//da qua inverto segno perchè i vincoli sono di maggiore uguale 
			for (GRBConstr c : constrDom) {
			    if (c.get(GRB.DoubleAttr.Slack) < 0.0) {
			         System.out.printf("[<" + c.get(GRB.StringAttr.ConstrName) + ">]");
			    }
			}
			for(GRBConstr[] riga : constrRag)
			    for (GRBConstr c : riga) 
				    if (c.get(GRB.DoubleAttr.Slack) < 0.0) 
				    	  System.out.printf("[<" + c.get(GRB.StringAttr.ConstrName) + ">]");
			System.out.printf("\n");	    
				      } 
		
		 /*
		  * per ogni vincolo del problema ottiene l'attributo PI e lo stampa a videp
		  */
		 private static void calcolaVariabiliDuale(GRBModel model)throws GRBException {
			int i=0;
			System.out.printf("la soluzione del duale è:\n");
			for (GRBConstr c : model.getConstrs()) {
			        System.out.printf("<lambda_%d>=<%.4f>, ",i,roundToFourthDecimal( c.get(GRB.DoubleAttr.Pi)));
			         i++;
			    
			}
		}
		
		 /*
		  * prende in input un valore h indice del magazzino
		  * e stampa a video gli attributi SARHSLow e SARHSUp relativi al vincolo
		  * di produzione collegato a quel magazzino
		  */
        private static void calcolaVariazioneSh(GRBModel model,int  h)throws GRBException{
			String string ="vincolo_produzione_i_"+h;
			GRBConstr c =model.getConstrByName(string);
			System.out.println("");
			System.out.println("intervallo s_"+h + "=["+roundToFourthDecimal(c.get(GRB.DoubleAttr.SARHSLow))+","+roundToFourthDecimal(c.get(GRB.DoubleAttr.SARHSUp))+"]");
		}
        
        /*
         * prende in input due valori x e y indici della variabile X e 
         * ne ottiene gli attributi SAObjLow e SAObjUp e divide entrambi per il prezzo unitario c
         * cosi da ottenere gli estremi di variazione della distanza tra il magazzino x e il cliente y
         * 
         */
		private static void calcolaVariazioneDxy(GRBModel model,int  x, int y )throws GRBException{
			GRBVar var=  model.getVarByName("xij_"+x+"_"+y);
			double d=roundToFourthDecimal(var.get(GRB.DoubleAttr.SAObjLow)/c);
			double d1=roundToFourthDecimal(var.get(GRB.DoubleAttr.SAObjUp)/c);
			if(d1>100000000000.0)
				System.out.println("intervallo d_"+x+"_"+y+" =["+d+"," +"INF)");	
			else System.out.println("intervallo d_"+x+"_"+y+" =["+d+"," +d1+"]");
		}
//end region
		
//start region TERZA PARTE==================================================================================================================
		private static void rispondiTerzaParte(GRBModel model, GRBVar []yi,GRBVar[] []xij,GRBModel rilassato)throws GRBException{
			System.out.println("\nQUESITO III:");
			
			model.set(GRB.IntParam.OutputFlag, 0);
			model.optimize();
			System.out.println("risparmio=<" + model.get(GRB.DoubleAttr.ObjVal)+">");
			
			stampaMagazzini(yi);
			magazzinoMenoUsato(xij);
			rilassato=model.relax();
			rilassato.set(GRB.IntParam.OutputFlag, 0);
			rilassato.optimize();
			System.out.println("rilassato continuo=<" + roundToFourthDecimal(rilassato.get(GRB.DoubleAttr.ObjVal))+">");
			
		}
		
		/*
		 * controlla i valori delle variabili binarie del problema di MPLI e stampa  gli indici
		 * dei magazzini collegate ad esse
		 * 
		 */
		private static void stampaMagazzini(GRBVar []yi)throws GRBException{

			System.out.printf("lista magazzini chiusi=[");
	     for (int i=0;i<yi.length; i++ )
	    	 if (yi[i].get(DoubleAttr.X )<0.98)
	    		 System.out.printf(" %d,",i);
	     System.out.printf(" ]");
					
	
} 
		
		/*
		 * prende in input le variabili xij del problema e ne estrae i valori
		 * calcola poi la somma della produzione per ciascun magazzino e trova poi i
		 * magazzini con la somma minore scartando quelli chiusi
		 * 
		 */
		public static void magazzinoMenoUsato(GRBVar[][]xij)throws GRBException{
			double[][]matrix=new double[xij.length][xij[0].length];
			for (int i = 0; i < xij.length; i++)
				for(int j=0; j<distanze[i].length; j++)
					matrix[i][j]=xij[i][j].get(DoubleAttr.X);
			ArrayList<Integer> minSumRows = findRowsWithMinSum(matrix);
		    System.out.printf("\nlista magazzini meno sfruttati=[");
		    for(int i: minSumRows)
		    	System.out.printf(" %d,", i);
		    System.out.println("]");
		}
		/*
		 * prende in input una matrice di double e ne fa la somma di ciascuna riga
		 * di quelle maggiori di zero ne calcola il minimo
		 */
		public static ArrayList<Integer> findRowsWithMinSum(double[][] matrix) {
		    // Inizializzazione dell'ArrayList che conterrà gli indici delle righe con la somma minore
		    ArrayList<Integer> minSumRows = new ArrayList<Integer>();

		    // Calcolo della somma di ciascuna riga e ricerca della somma minore
		    double minSum = Double.MAX_VALUE;
		    for (int i = 0; i < matrix.length; i++) {
		        double rowSum = 0;
		        for (int j = 0; j < matrix[i].length; j++) {
		            rowSum += matrix[i][j];
		        }
		       if(rowSum>0.0001) {
		        if (Math.abs(rowSum - minSum)<0.0001 && (rowSum>0.0001)) {
		           minSumRows.add(i);
		        } else if(rowSum <minSum) {
		        	  minSumRows.clear();
			            minSumRows.add(i);
			            minSum = rowSum;
		        }
		       }
		    }

		    // Restituzione dell'ArrayList degli indici delle righe con la somma minore
		    return minSumRows;
		}

		//CREAZIONE MODELLO
    	private static void impostaParametri(GRBEnv env) throws GRBException 
		{
			env.set(GRB.IntParam.Method, 0);
			env.set(GRB.IntParam.Presolve, 0);
		}
		private static GRBVar[][] aggiungiVariabili(GRBModel model, int[] produzione, int[] domanda) throws GRBException 
		{
			GRBVar[][] xij = new GRBVar[produzione.length][domanda.length];
			
			for (int i = 0; i < produzione.length; i++) 
			{
				for (int j = 0; j < domanda.length; j++) 
				{
					xij[i][j] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "xij_"+i+"_"+j);
					
				}
			}
			return xij;
		}
		private static void aggiungiFunzioneObiettivo(GRBModel model, GRBVar[][] xij, int[][] distanze, double c) throws GRBException 
		{
			GRBLinExpr obj = new GRBLinExpr();
			
			for (int i = 0; i < distanze.length; i++) 
			{
				for (int j = 0; j < distanze[i].length; j++) 
				{
					obj.addTerm(distanze[i][j]*c, xij[i][j]);
				}
			}
			model.setObjective(obj);
			model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);
		}
	
		private static GRBConstr[] aggiungiVincoliProduzione(GRBModel model, GRBVar[][] xij, int[] produzione) throws GRBException 
		{
			GRBConstr[] xi = new GRBConstr[produzione.length];
			for (int i = 0; i < produzione.length; i++) 
			{
				GRBLinExpr expr = new GRBLinExpr();
				
				for (int j = 0; j < xij[0].length; j++) 
				{
					expr.addTerm(1, xij[i][j]);
				}
				xi[i] = model.addConstr(expr, GRB.LESS_EQUAL, produzione[i], "vincolo_produzione_i_"+i);
			}
			return xi;
		}
		private static GRBConstr[]aggiungiVincoliDomanda(GRBModel model, GRBVar[][] xij, int[] domanda) throws GRBException 
		{
			GRBConstr[] xj = new GRBConstr[domanda_clienti.length];
			for (int j = 0; j < domanda.length; j++) 
			{
				GRBLinExpr expr = new GRBLinExpr();
				
				for (int i = 0; i < xij.length; i++) 
				{
					expr.addTerm(1, xij[i][j]);
				}
				xj[j] = model.addConstr(expr, GRB.GREATER_EQUAL, domanda[j], "vincolo_domanda_j_"+j);
			}
			return xj;
		}
		private static GRBConstr[][] aggiungiVincoliRaggio(GRBModel model, GRBVar[][] xij, int[][] distanze,int  k) throws GRBException 
		{
			GRBConstr[][] constrij = new GRBConstr[produzione.length][domanda_clienti.length];
			for (int i = 0; i < distanze.length; i++)
				for(int j=0; j<distanze[i].length; j++)
			{
				GRBLinExpr expr = new GRBLinExpr();
				expr.addTerm((k-distanze[i][j]), xij[i][j]);
				constrij[i][j]=model.addConstr(expr, GRB.GREATER_EQUAL, 0, "vincolo_raggio_ij_"+i+"_"+j);
			}
			return constrij;
		}
		private static void modificaVincoliRaggio(GRBModel model, GRBVar[][] xij, int[][] distanze,int  k) throws GRBException 
		{
			for (int i = 0; i < distanze.length; i++)
				for(int j=0; j<distanze[i].length; j++) {
				GRBConstr c1 = model.getConstrByName("vincolo_raggio_ij_"+i+"_"+j);
				model.remove(c1);
				GRBLinExpr expr = new GRBLinExpr();
				expr.addTerm((k-distanze[i][j]), xij[i][j]);
				model.addConstr(expr, GRB.GREATER_EQUAL, 0, "vincolo_raggio_ij_"+i+"_"+j);
				
			}
		}
	
		
		private static void aggiungiVincoliProduzioneIntera(GRBModel modelIntero, GRBVar[][] xij, GRBVar[] yi,
				int[] produzione) throws GRBException{
			for (int i = 0; i < produzione.length; i++) {
			
				GRBLinExpr exprD = new GRBLinExpr();
				GRBLinExpr exprS = new GRBLinExpr();
				
				
				for (int j = 0; j < xij[0].length; j++) 
				{
					exprS.addTerm(1, xij[i][j]);
				}
				
				   exprD.addTerm(produzione[i], yi[i]);
				
				modelIntero.addConstr(exprS, GRB.LESS_EQUAL, exprD, "vincolo_produzione_i_"+i);
			}
			
		}
		private static void aggiungiFunzioneObiettivoIntera(GRBModel modelIntero, GRBVar[] yi, int[] alfa_j) throws GRBException{
               GRBLinExpr obj = new GRBLinExpr();
			
			for (int i = 0; i < alfa_j.length; i++) {
				obj.addConstant(alfa_j[i]);
			    obj.addTerm(-alfa_j[i], yi[i]);
				}
			    modelIntero.setObjective(obj, GRB.MAXIMIZE);
			
			
		}
        private static GRBVar[] aggiungiVariabiliIntere(GRBModel modelIntero, int[] alfa_j) throws GRBException{
           GRBVar[] alfaj = new GRBVar[alfa_j.length];
			
			for (int i = 0; i < alfa_j.length; i++) {
				
			alfaj[i] = modelIntero.addVar(0,  GRB.INFINITY,0,GRB.BINARY, "yi_"+i);
			}
			return alfaj;
         
		}

		
        //LETTURA FILE	
        public static double roundToFourthDecimal(double num) {
            double roundedNum = Math.round(num * 10000.0) / 10000.0;
            return roundedNum;
        }
			     public static void fileToVariable() {
			         String fileName = "singolo_13.txt";

			         try {
			             processFile(fileName);
			         } catch (IOException e) {
			             e.printStackTrace();
			         }
			     }
			     public static void processFile(String fileName) throws IOException {
			    	 
			     boolean leggendor_i=false;
			     boolean leggendos_j=false;
			     boolean leggendoalfa_j=false;
			     int i=0;
			     
			         try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
			             String line;
			             

			             while ((line = br.readLine()) != null) {
			            	 String[] intElements = line.split(" ");
			             
			                 // Esempio: leggere un array di interi
			            	
			            	if (leggendoMatricedij) 
			            	{
			            		for (int j = 0; j < n; j++) {
			                        distanze[i][j] = Integer.parseInt(intElements[j]);
			                    }
			                    i++;
			                    if(i>=m)
			                    	leggendoMatricedij=false;
			            	}
			            	if (leggendos_j) 
			            	{
			            		for (int j = 0; j < m; j++) {
			                        produzione[j] = Integer.parseInt(intElements[j]);
			                    }
			                    
			                    	leggendos_j=false;
			            	}
			            	if (leggendor_i) 
			            	{
			            		for (int j = 0; j < n; j++) {
			                        domanda_clienti[j] = Integer.parseInt(intElements[j]);
			                    }
			                    
			                    	leggendor_i=false;
			            	}
			            	if (leggendoalfa_j) 
			            	{
			            		for (int j = 0; j < m; j++) {
			                        alfa_j[j] = Integer.parseInt(intElements[j]);
			                    }
			                    
			                    	leggendoalfa_j=false;
			            	}
			            	 if (intElements[0].equals("n")) {
			            		 n=(Integer.parseInt(intElements[1]));
			            		 
			            	 }
			            	 if (intElements[0].equals("m")) {
			            		 m=(Integer.parseInt(intElements[1]));
			            	 }
			            	 if (intElements[0].equals("c")) {
			            		 c=(Double.parseDouble(intElements[1]));
			            	 }
			            	 if (intElements[0].equals("k")) {
			            		 k=(Integer.parseInt(intElements[1]));
			            	 }
			            	 if (intElements[0].equals("h")) {
			            		 h=(Integer.parseInt(intElements[1]));
			            	 }
			            	 if (intElements[0].equals("x")) {
			            		 x=(Integer.parseInt(intElements[2]));
			            	 }
			            	 if (intElements[1].equals("y")) {
			            		 y=(Integer.parseInt(intElements[3]));
			            	 }
			            	 if (intElements[0].equals("matrice")) {
			            		distanze=new int [m][n];
			            		leggendoMatricedij=true;
			            		i=0;
			            	 }
			            	 if (intElements[0].equals("r_i")) {
			            		    domanda_clienti= new int [n];
			            		    leggendor_i=true;
				            		i=0;
				            	 }
			            	 if (intElements[0].equals("s_j")) {
			            		    produzione= new int [m];
				            		leggendos_j=true;
				            		i=0;
				            	 }
			            	 if (intElements[0].equals("alfa_j")) {
			            		    alfa_j= new int [m];
				            		leggendoalfa_j=true;
				            		i=0;
				            	 }
			            	
			            		 
			             }
			                    
			         }
			     }
			     private static void stampaVettore(int[] vettore) {
	    	  for(int i = 0; i < vettore.length; i++) {
	    	        System.out.printf( "%d   ", vettore[i]);
	    	    }
	    }
			     private static void stampaVettore(double[] vettore) {
			    	  for(int i = 0; i < vettore.length; i++) {
			    	        System.out.printf( "%.2f  ", vettore[i]);
			    	    }
			    }
			     private static void stampaMatrice(int[][] vettore) {
	    	 
	   	  for(int i = 0; i < vettore.length; i++) {
	   		System.out.printf("\n ");
	   		  for(int j=0; j<vettore[i].length;j++)
	   	        System.out.printf("%d ", vettore[i][j]);
	   	    }
	    }
	}