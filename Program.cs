using System;

namespace Projekt1
{
	static class Globals
	{
		public static double H = 0.100; //wysokosc
		public static double L = 0.100; //dlugosc
		public static int nH = 4; //ilosc elementow w wysokosci
		public static int nL = 4; //ilosc elementow w dlugosci
		public static int nW = nH * nL; //ilosc wezlow
		public static int E = (nH - 1) * (nL - 1); //ilosc elementow

		public static int weight = 1;
		public static double conductivity = 25;
		public static double convection = 300;
		public static double ambientTemp = 1200; //temperatura otoczenia
		public static double specificHeat = 700; //cieplo wlasciwe
		public static double density = 7800; //gestosc
		public static double initialTemperature = 100;
		public static int simulationTime = 500;
		public static int timeStep = 25;
		public static double plus = (1 / Math.Sqrt(3));
		public static double minus = -(1 / Math.Sqrt(3));

		public static double EPSILON = 1e-10;

	}

	static class GaussianElimination
	{

		public static double[]  lsolve(double[,] A2, double[] b2)
		{

			int n = Globals.nW;
			double[,] A = new double[n, n];
			double[] b = new double[n];


			for (int i = 0; i < n; i++)
			{
				b[i] = b2[i];
				for (int j = 0; j < n; j++)
				{
					A[i, j] = A2[i, j];
				}
			}


			for (int p = 0; p < n; p++)
			{

				int max = p;
				for (int i = p + 1; i < n; i++)
				{
					if (Math.Abs(A[i, p]) > Math.Abs(A[max, p]))
					{
						max = i;
					}
				}


				double[] temp = new double[n];
				for (int d = 0; d < n; d++)
				{
					temp[d] = A[p, d];
					A[p,d] = A[max,d];
					A[max,d] = temp[d];

					double t = b[p];
					b[p] = b[max];
					b[max] = t;
				}
				
				if (Math.Abs(A[p, p]) <= Globals.EPSILON)
				{
					return null;
				}

				for (int i = p + 1; i < n; i++)
				{
					double alpha = A[i,p] / A[p,p];
					b[i] -= alpha * b[p];
					for (int j = p; j < n; j++)
					{
						A[i,j] -= alpha * A[p,j];
					}
				}
			}

			double[] x = new double[n];
			for (int i = n - 1; i >= 0; i--)
			{
				double sum = 0.0;
				for (int j = i + 1; j < n; j++)
				{
					sum += A[i,j] * x[j];
				}
				x[i] = (b[i] - sum) / A[i,i];
			}

		return x;
		}
	}



	class Node
	{
		public double x = 0;
		public double y = 0;
		public double t = 0;
		public int bc = 0;
		public int id = 0;
		public void show()
		{
			Console.Write("X: " + x + " Y: " + y + " BC: " + bc);
		}

	};

	class Element
	{
		public Node[] points;

		public double[,] derTab = new double[4, 4]; 
		public double[] jakobians = new double[4]; 
		public double[,] dNdx = new double[4, 4];
		public double[,] dNdy = new double[4, 4];

		public double[,,] matrixHdx1 = new double[4, 4, 4]; 
		public double[,,] matrixHdy1 = new double[4, 4, 4]; 
		public double[,,] matrixH2 = new double[4, 4, 4]; 
		public double[,] matrixH = new double[4, 4]; 

		public Surface[] surface = new Surface[4]; 
		public double[,] matrixH_BC = new double[4, 4];
		public double[] vectorP = new double[4];
		public double[,] matrixC = new double[4, 4];


		public void show()
		{
			//Console.Write("Gorny lewy");
			points[3].show();
			Console.Write(" \t\t\t   ||  ");
			//Console.Write("Gorny prawy");
			points[2].show();

			Console.WriteLine();

			//Console.Write("Dolny Lewy ");
			points[0].show();
			//Console.Write(" , Dolny Prawy ");
			Console.Write(" \t\t\t   ||  ");
			points[1].show();
			Console.WriteLine();


		}
	};

	//--------------------------------------- POWIERZCHNIA (BOK ELEMENTU) ------------------------------------------------------
	public class Surface
	{
		public int bc;
		public double[] ksi = new double[2];
		public double[] eta = new double[2];
		public double[,] N = new double[2, 4]; 
		public double detJ;
		public double[,] partForH = new double[4, 4];
		public double[] partForP = new double[4];
	};

	// --------------------------------------- SIATKA -----------------------------------------------------------------------------------
	class siatka
	{
		Element[] elements;
		Node[] nodes;
		double[,] globalH;
		double[,] globalC;
		double[] globalP;
		double[] currentGlobalP; 

		double[] temperatures;

		public siatka(int element, int node)
		{
			elements = new Element[element + 1];
			nodes = new Node[node + 1];

			globalP = new double[Globals.nW];
			currentGlobalP = new double[Globals.nW];
			globalH = new double[Globals.nW, Globals.nW];
			globalC = new double[Globals.nW, Globals.nW];

		}


		public void fillNodes()
		{
			double dX = Globals.L / (Globals.nL - 1);
			double dY = Globals.H / (Globals.nH - 1);

			int currentNode = 1;
			for (int i = 1; i <= Globals.nL; i++)
			{
				for (int j = 1; j <= Globals.nH; j++)
				{
					nodes[currentNode] = new Node();
					nodes[currentNode].x = dX * (i - 1);
					nodes[currentNode].y = dY * (j - 1);


					nodes[currentNode].id = currentNode - 1;
					Console.WriteLine("NODE: " + currentNode + " X: " + nodes[currentNode].x + " Y: " + nodes[currentNode].y);

					currentNode++;
				}
			}

			
			for (int i = 1; i <= Globals.nW; i++)
			{
				if ((i >= 1 && i <= Globals.nH) || (i > Globals.nW - Globals.nH && i <= Globals.nW) || ((i - 1) % Globals.nH == 0) || i % Globals.nH == 0)
				{
					nodes[i].bc = 1;
				}
				else
				{
					nodes[i].bc = 0;
				}
			}
		}
		public void fillElements()
		{

			int fix = 0;
			int fixCount = 0;
			for (int i = 1; i <= Globals.E; i++)

			{
				elements[i] = new Element();
				elements[i].points = new Node[4];
				elements[i].points[0] = nodes[i + fix];
				elements[i].points[1] = nodes[i + Globals.nH + fix];
				elements[i].points[2] = nodes[i + Globals.nH + 1 + fix];
				elements[i].points[3] = nodes[i + 1 + fix];
				Console.WriteLine("ELEMENT " + i + " ");
				elements[i].show();


				fixCount++;
				if (fixCount == Globals.nH - 1)
				{
					fixCount = 0;
					fix += 1;
				}
			}
		}

		public void calculateDeriatives(UniversalEl uniEl)
		{

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						elements[currentElement].derTab[i, j] = 0;
					}
				}


				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{
					for (int i = 0; i < 4; i++)
					{
						//dxde
						elements[currentElement].derTab[currentPoint, 0] += uniEl.tabDerKsiN[currentPoint, i] * elements[currentElement].points[i].x;
						//dyde
						elements[currentElement].derTab[currentPoint, 1] += uniEl.tabDerKsiN[currentPoint, i] * elements[currentElement].points[i].y;
						//dxdn
						elements[currentElement].derTab[currentPoint, 2] += uniEl.tabDerEtaN[currentPoint, i] * elements[currentElement].points[i].x;
						//dydn
						elements[currentElement].derTab[currentPoint, 3] += uniEl.tabDerEtaN[currentPoint, i] * elements[currentElement].points[i].y;
					}
				}

				Console.WriteLine("ELEMENT:" + currentElement);
				for (int i = 0; i < 4; i++)
				{
					Console.Write("POINT:" + i + ": ");
					for (int j = 0; j < 4; j++)
					{
						switch (j)
						{
							case 0:
								Console.Write("dxdE:" + elements[currentElement].derTab[i, j] + " ");
								break;
							case 1:
								Console.Write("dydE:" + elements[currentElement].derTab[i, j] + " ");
								break;
							case 2:
								Console.Write("dxdn:" + elements[currentElement].derTab[i, j] + " ");
								break;
							case 3:
								Console.Write("dydn:" + elements[currentElement].derTab[i, j] + " ");
								break;
						}
					}
					Console.WriteLine();
				}
				Console.WriteLine();
			}
		}

		public void calculateJakobians()
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{
					double a = elements[currentElement].derTab[currentPoint, 0];
					double b = elements[currentElement].derTab[currentPoint, 1];
					double c = elements[currentElement].derTab[currentPoint, 2];
					double d = elements[currentElement].derTab[currentPoint, 3];
					elements[currentElement].jakobians[currentPoint] = (a * d - b * c);
				}

			}
		}

		public void calculateDeriatives2(UniversalEl uniEl)
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				double[,] tab = new double[4, 4];

				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{

					tab[currentPoint, 0] = (elements[currentElement].derTab[currentPoint, 3]) / elements[currentElement].jakobians[currentPoint];
					tab[currentPoint, 1] = (-elements[currentElement].derTab[currentPoint, 1]) / elements[currentElement].jakobians[currentPoint];
					tab[currentPoint, 2] = (-elements[currentElement].derTab[currentPoint, 2]) / elements[currentElement].jakobians[currentPoint];
					tab[currentPoint, 3] = (elements[currentElement].derTab[currentPoint, 0]) / elements[currentElement].jakobians[currentPoint];
				}
			
				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{
					for (int i = 0; i < 4; i++)
					{
						elements[currentElement].dNdx[currentPoint, i] = tab[currentPoint, 0] * uniEl.tabDerKsiN[currentPoint, i] + tab[currentPoint, 1] * uniEl.tabDerEtaN[currentPoint, i];
						elements[currentElement].dNdy[currentPoint, i] = tab[currentPoint, 2] * uniEl.tabDerKsiN[currentPoint, i] + tab[currentPoint, 3] * uniEl.tabDerEtaN[currentPoint, i];
					}
				}

			}
		}

		public void H_step1()
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{

				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{

					for (int i = 0; i < 4; i++)
					{
						for (int j = 0; j < 4; j++)
						{
							elements[currentElement].matrixHdx1[currentPoint, i, j] = elements[currentElement].dNdx[currentPoint, i] * elements[currentElement].dNdx[currentPoint, j];
							elements[currentElement].matrixHdy1[currentPoint, i, j] = elements[currentElement].dNdy[currentPoint, i] * elements[currentElement].dNdy[currentPoint, j];
						}
					}


				}
			}
		}

		public void H_step2()
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{


				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{

					for (int i = 0; i < 4; i++)
					{
						for (int j = 0; j < 4; j++)
						{
							elements[currentElement].matrixHdx1[currentPoint, i, j] *= elements[currentElement].jakobians[currentPoint];
							elements[currentElement].matrixHdy1[currentPoint, i, j] *= elements[currentElement].jakobians[currentPoint];
						}
					}


				}
			}
		}

		public void H_step3()
		{
			double k = Globals.conductivity;

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int currentPoint = 0; currentPoint < 4; currentPoint++)
				{

					for (int i = 0; i < 4; i++)
					{
						for (int j = 0; j < 4; j++)
						{
							elements[currentElement].matrixH2[currentPoint, i, j] = k * (elements[currentElement].matrixHdx1[currentPoint, i, j] + elements[currentElement].matrixHdy1[currentPoint, i, j]);
						}
					}

				}
			}
		}

		public void H_with_no_bc()
		{

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				Console.WriteLine("\n ELEMENT: " + currentElement);
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						elements[currentElement].matrixH[i, j] = elements[currentElement].matrixH2[0, i, j] + elements[currentElement].matrixH2[1, i, j] + elements[currentElement].matrixH2[2, i, j] + elements[currentElement].matrixH2[3, i, j];
					}
				}

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						Console.Write(elements[currentElement].matrixH[i, j] + " ");
					}
					Console.WriteLine();


				}
			}
		}

		public void calculate_bc()
		{
			
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int i = 0; i < 4; i++) {
					elements[currentElement].surface[i] = new Surface();
					elements[currentElement].surface[i].bc = 0;
				} 

				if (elements[currentElement].points[0].bc == 1 && elements[currentElement].points[1].bc == 1) elements[currentElement].surface[0].bc = 1;
				if (elements[currentElement].points[1].bc == 1 && elements[currentElement].points[2].bc == 1) elements[currentElement].surface[1].bc = 1;
				if (elements[currentElement].points[2].bc == 1 && elements[currentElement].points[3].bc == 1) elements[currentElement].surface[2].bc = 1;
				if (elements[currentElement].points[3].bc == 1 && elements[currentElement].points[0].bc == 1) elements[currentElement].surface[3].bc = 1;


				
			}
		}

		public void setLocalCoordsOnSurfaces()
		{
			
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				
				elements[currentElement].surface[0].ksi[0] = Globals.minus;
				elements[currentElement].surface[0].ksi[1] = Globals.plus;
				elements[currentElement].surface[0].eta[0] = -1;
				elements[currentElement].surface[0].eta[1] = -1;
			
				elements[currentElement].surface[1].ksi[0] = 1;
				elements[currentElement].surface[1].ksi[1] = 1;
				elements[currentElement].surface[1].eta[0] = Globals.minus;
				elements[currentElement].surface[1].eta[1] = Globals.plus;

				elements[currentElement].surface[2].ksi[0] = Globals.plus;
				elements[currentElement].surface[2].ksi[1] = Globals.minus;
				elements[currentElement].surface[2].eta[0] = 1;
				elements[currentElement].surface[2].eta[1] = 1;

				elements[currentElement].surface[3].ksi[0] = -1;
				elements[currentElement].surface[3].ksi[1] = -1;
				elements[currentElement].surface[3].eta[0] = Globals.plus;
				elements[currentElement].surface[3].eta[1] = Globals.minus;
			}
		}

		public void calculate_N_onSurfaces()
		{

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int currentSurface = 0; currentSurface < 4; currentSurface++)
				{
					for (int currentPoint = 0; currentPoint < 2; currentPoint++)
					{
						for (int currentFunction = 0; currentFunction < 4; currentFunction++)
						{
							elements[currentElement].surface[currentSurface].N[currentPoint, currentFunction] = Function_of_Shape.N(
								elements[currentElement].surface[currentSurface].ksi[currentPoint],
								elements[currentElement].surface[currentSurface].eta[currentPoint],
								currentFunction);
						}
					}

				}
			}

		}

		public double jakobian(double x1, double x2, double y1, double y2)
		{
		
			return (Math.Sqrt(Math.Pow(x2 - x1, 2) + Math.Pow(y2 - y1, 2))) / 2;
		}

		public void calculateJakobiansOnSurfaces()
		{

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{ 
				elements[currentElement].surface[0].detJ = jakobian(elements[currentElement].points[0].x, elements[currentElement].points[1].x, elements[currentElement].points[0].y, elements[currentElement].points[1].y);
				elements[currentElement].surface[1].detJ = jakobian(elements[currentElement].points[1].x, elements[currentElement].points[2].x, elements[currentElement].points[1].y, elements[currentElement].points[2].y);
				elements[currentElement].surface[2].detJ = jakobian(elements[currentElement].points[2].x, elements[currentElement].points[3].x, elements[currentElement].points[2].y, elements[currentElement].points[3].y);
				elements[currentElement].surface[3].detJ = jakobian(elements[currentElement].points[3].x, elements[currentElement].points[0].x, elements[currentElement].points[3].y, elements[currentElement].points[0].y);

			}
		}

		public void calculateH_BC_1()
		{
			
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{

				for (int currentSurface = 0; currentSurface < 4; currentSurface++)
				{

					for (int i = 0; i < 4; i++)
					{
						for (int j = 0; j < 4; j++)
						{
							double pc1part = elements[currentElement].surface[currentSurface].N[0, i] * elements[currentElement].surface[currentSurface].N[0, j] * Globals.convection;
							double pc2part = elements[currentElement].surface[currentSurface].N[1, i] * elements[currentElement].surface[currentSurface].N[1, j] * Globals.convection;

							elements[currentElement].surface[currentSurface].partForH[i, j] = (pc1part + pc2part) * elements[currentElement].surface[currentSurface].detJ;
						}
						
					}
				}
			}
		}

		public void calculateH_BC_2()
		{

			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{

				
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						elements[currentElement].matrixH_BC[i,j] =
																	elements[currentElement].surface[0].bc * elements[currentElement].surface[0].partForH[i,j] +
																	elements[currentElement].surface[1].bc * elements[currentElement].surface[1].partForH[i,j] +
																	elements[currentElement].surface[2].bc * elements[currentElement].surface[2].partForH[i,j] +
																	elements[currentElement].surface[3].bc * elements[currentElement].surface[3].partForH[i,j];
						
					}
					
				}
			}
		}

		public void calculate_H_final()
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						elements[currentElement].matrixH[i,j] += elements[currentElement].matrixH_BC[i,j];
					
					}
					
				}
			}
		}

		public void P()
		{
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				for (int currentSurface = 0; currentSurface < 4; currentSurface++)
				{
					for (int i = 0; i < 4; i++)
					{
						double pc1part = elements[currentElement].surface[currentSurface].N[0,i] * Globals.convection * Globals.ambientTemp;
						double pc2part = elements[currentElement].surface[currentSurface].N[1,i] * Globals.convection * Globals.ambientTemp;
						elements[currentElement].surface[currentSurface].partForP[i] = (pc1part + pc2part) * elements[currentElement].surface[currentSurface].detJ;
					}
				} 

				for (int i = 0; i < 4; i++)
				{
					elements[currentElement].vectorP[i] = elements[currentElement].surface[0].partForP[i] * elements[currentElement].surface[0].bc +
														  elements[currentElement].surface[1].partForP[i] * elements[currentElement].surface[1].bc +
														  elements[currentElement].surface[2].partForP[i] * elements[currentElement].surface[2].bc +
														  elements[currentElement].surface[3].partForP[i] * elements[currentElement].surface[3].bc;
				}

			
			}
		}

		public void C(UniversalEl uniEl)
		{
			
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{
				
				double[,,]partialC=new double [4,4,4];  
				for (int currentPoint = 0; currentPoint< 4; currentPoint++) {
					
					for (int i = 0; i< 4; i++) {
						for (int j = 0; j< 4; j++) {
							partialC[currentPoint,i,j] = uniEl.tabN[currentPoint,i]* uniEl.tabN[currentPoint,j] * elements[currentElement].jakobians[currentPoint] * Globals.specificHeat* Globals.density;
						}
					}
				}
		
				for (int i = 0; i< 4;i++) {
					for (int j = 0; j< 4;j++) {
						elements[currentElement].matrixC[i,j] = partialC[0,i,j] + partialC[1,i,j] + partialC[2,i,j] + partialC[3,i,j];
					}
					
				}

			}
		}

		public void agregate()
		{
		
			for (int i = 0; i < Globals.nW; i++)
			{
				globalP[i] = 0;
				for (int j = 0; j < Globals.nW; j++)
				{
					globalH[i,j] = 0;
					globalC[i,j] = 0;
				}
			}
			for (int currentElement = 1; currentElement < Globals.E + 1; currentElement++)
			{

				int[] ID= new int[4]; 

				for (int j = 0; j < 4; j++)
				{
					ID[j] = elements[currentElement].points[j].id; 
				}



				for (int k = 0; k < 4; k++)
				{
					globalP[ID[k]] += elements[currentElement].vectorP[k];
					for (int z = 0; z < 4; z++)
					{
						globalH[ID[k],ID[z]] += elements[currentElement].matrixH[k,z];
						globalC[ID[k],ID[z]] += elements[currentElement].matrixC[k,z];
					}

				}

			}

			
			for (int i = 0; i < Globals.nW; i++)
			{
				for (int j = 0; j < Globals.nW; j++)
				{
					Console.WriteLine ( globalH[i,j] + " ");
				}
					Console.WriteLine ();
			}

			//GLOBAL C:";
			for (int i = 0; i < Globals.nW; i++)
			{
				for (int j = 0; j < Globals.nW; j++)
				{
					Console.Write(globalC[i,j] +" ");
				}
				Console.WriteLine();
			}

			//LOBAL P: 
			for (int i = 0; i < Globals.nW; i++)
			{
				Console.Write( globalP[i] + " ");
			}
			Console.WriteLine();
		}

		public void makeUltimateH()
		{

			for (int i = 0; i < Globals.nW; i++)
			{
				for (int j = 0; j < Globals.nW; j++)
				{
					globalC[i,j] = globalC[i,j] / Globals.timeStep;
				}
			}
			for (int i = 0; i < Globals.nW; i++)
			{
				for (int j = 0; j < Globals.nW; j++)
				{
					globalH[i,j] = globalH[i,j] + globalC[i,j];
					Console.Write ( globalH[i,j] + " ");
				}
				Console.WriteLine ();
			}
		}

		public void makeUltimateP(double[] T0)
		{
			int L = Globals.nW;
			double[] partialP = new double[L];
			for (int i = 0; i < L; i++) partialP[i] = 0;

			
			for (int i = 0; i < L; i++)
			{
				for (int j = 0; j < L; j++)
				{
					partialP[i] += (globalC[i,j] * T0[j]);
				}
			}

			
			for (int i = 0; i < L; i++)
			{
				currentGlobalP[i] = globalP[i] + partialP[i]; 		
			}
			
		}


		public void findTemperatures()
		{
			temperatures = GaussianElimination.lsolve(globalH, currentGlobalP);
			
		}

		public double findMin(double[] A)
		{

			double min = 9999999999;

			for (int i = 0; i < Globals.nW; i++)
			{
				if (A[i] < min) min = A[i];
			}
			return min;
		}

		double findMax(double[] A)
		{
			double max = -9999999999;

			for (int i = 0; i < Globals.nW; i++)
			{
				if (A[i] > max) max = A[i];
			}
			return max;
		}
		public void wyniki()
		{

			makeUltimateH();

			for (int i = Globals.timeStep; i <= Globals.simulationTime; i += Globals.timeStep)
			{
				
				if (i == Globals.timeStep)
				{ 

					int x = Globals.nW;
					double[] T0 = new double[x];
					for (int j = 0; j < x; j++) T0[j] = Globals.initialTemperature;

					makeUltimateP(T0); 

					findTemperatures(); 
					Console.WriteLine("Time[s]:" + i + " MinTemp:" + findMin(temperatures) + " MaxTemp:" + findMax(temperatures) );
					continue;
				}

				makeUltimateP(temperatures); 

				findTemperatures();
				Console.WriteLine("Time[s]:" + i + " MinTemp:" + findMin(temperatures) + " MaxTemp:" + findMax(temperatures));
			}
		}
	}
	//--------------------- FUNKCJE KSZTAŁTU -----------------------------------------------------------
	static class Function_of_Shape {
			public static double N(double ksi, double eta, int functionNumber)
			{
				switch (functionNumber)
				{
					case 0:
						return 0.25 * (1 - ksi) * (1 - eta);  //N1
					case 1:
						return 0.25 * (1 + ksi) * (1 - eta); //N2
					case 2:
						return 0.25 * (1 + ksi) * (1 + eta); //N3
					case 3:
						return 0.25 * (1 - ksi) * (1 + eta); //N4
					default:
						//nothing
						return 0;
				}

			}

			public static double derKsiN(double eta, int functionNumber)
			{  

				switch (functionNumber)
				{
					case 0:
						return -0.25 * (1 - eta);  //N1
					case 1:
						return 0.25 * (1 - eta); //N2
					case 2:
						return 0.25 * (1 + eta); //N3
					case 3:
						return -0.25 * (1 + eta);  //N4
					default:
						//nothing
						return 0;
				}
			}

			public static double derEtaN(double ksi, int functionNumber)
			{

				switch (functionNumber)
				{
					case 0:
						return -0.25 * (1 - ksi);  //N1
					case 1:
						return -0.25 * (1 + ksi); //N2
					case 2:
						return 0.25 * (1 + ksi); //N3
					case 3:
						return 0.25 * (1 - ksi);  //N4
					default:
						//nothing
						return 0;
				}
			}
		}
	
		//------------------- TWORZENIE ELEMENTU UNIWERSALNEGO -------------------------------------------------

		class Point
		{
			public double ksi;
			public double eta;

			public Point() {
				ksi = 0;
				eta = 0;
			}

			public Point(double ksi_lok, double eta_lok) {
				ksi = ksi_lok;
				eta = eta_lok;
			}

		};

		class UniversalEl
		{
			public int weight1;
			public int weight2;
			public int weight3;
			public int weight4;
			public Point[] points = new Point[4];

			public double[,] tabN;
			public double[,] tabDerKsiN;
			public double[,] tabDerEtaN;
			public UniversalEl() {
			}
			public UniversalEl(int weight1_l, int weight2_l, int weight3_l, int weight4_l, Point one, Point two, Point three, Point four) {
				weight1 = weight1_l;
				weight2 = weight2_l;
				weight3 = weight3_l;
				weight4 = weight4_l;

				points[0] = new Point();
				points[0] = one;
				points[1] = new Point();
				points[1] = two;
				points[2] = new Point();
				points[2] = three;
				points[3] = new Point();
				points[3] = four;


				tabN = createTableN();
				tabDerKsiN = createTableDerKsiN();
				tabDerEtaN = createTableDerEtaN();
			}
			public void show()
			{
				Console.WriteLine("Printing Universal Element 2 x 2:");
				Console.WriteLine("Point 1:  ksi: " + points[0].ksi + " eta: " + points[0].eta);
				Console.WriteLine("Point 2:  ksi: " + points[1].ksi + " eta: " + points[1].eta);
				Console.WriteLine("Point 3:  ksi: " + points[2].ksi + " eta: " + points[2].eta);
				Console.WriteLine("Point 4:  ksi: " + points[3].ksi + " eta: " + points[3].eta);

			}
			public UniversalEl createUniEl2x2()
			{

				Point one = new Point(Globals.minus, Globals.minus);
				Point two = new Point(Globals.plus, Globals.minus);
				Point three = new Point(Globals.plus, Globals.plus);
				Point four = new Point(Globals.minus, Globals.plus);
				UniversalEl a = new UniversalEl(Globals.weight, Globals.weight, Globals.weight, Globals.weight, one, two, three, four);

				return a;
			}


			public double[,] createTableN()
			{
				double[,] tabN;
				tabN = new double[4, 4];

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						Point currentPoint = points[i];
						 
						double result = Function_of_Shape.N(currentPoint.ksi, currentPoint.eta, j);
						tabN[i, j] = result;
					}
				}

				return tabN;
			}

			public double[,] createTableDerKsiN()
			{

				double[,] tabN;
				tabN = new double[4, 4];

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						Point currentPoint = points[i];
						
						double result = Function_of_Shape.derKsiN(currentPoint.eta, j);
						tabN[i, j] = result;
					}
				}

				return tabN;
			}

			public double[,] createTableDerEtaN()
			{

				double[,] tabN;
				tabN = new double[4, 4];

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						Point currentPoint = points[i];
						
						double result = Function_of_Shape.derEtaN(currentPoint.ksi, j);
						tabN[i, j] = result;
					}
				}

				return tabN;
			}
		}

	class Program
		{
			static void printTab4x4(double[,] tab)
			{
				for (int i = 0; i < 4; i++)
				{
					Console.Write("PC " + i + ": ");
					for (int j = 0; j < 4; j++)
					{
						Console.Write(tab[i, j] + "   ");
					}
					Console.WriteLine();
				}
			}
			static void Main(string[] args)
			{
				siatka grid = new siatka(Globals.E, Globals.nW);
				grid.fillNodes();
				grid.fillElements();

				Console.WriteLine();
				UniversalEl uni = new UniversalEl();
				uni = uni.createUniEl2x2();
				uni.show();

				Console.WriteLine();
				Console.WriteLine(" N:");

				printTab4x4(uni.tabN);
				Console.WriteLine();

				Console.WriteLine(" KSI OF N:");
				printTab4x4(uni.tabDerKsiN);

				Console.WriteLine(" ETA OF N:");
				printTab4x4(uni.tabDerEtaN);
				Console.WriteLine();

				grid.calculateDeriatives(uni);
				grid.calculateJakobians();
				grid.calculateDeriatives2(uni);

				grid.H_step1();
				grid.H_step2();
				grid.H_step3();
				grid.H_with_no_bc();

			
				grid.calculate_bc();
				grid.setLocalCoordsOnSurfaces();
				grid.calculate_N_onSurfaces();

				grid.calculateJakobiansOnSurfaces();
				grid.calculateH_BC_1();
				grid.calculateH_BC_2();
				grid.calculate_H_final();

				grid.P();
				grid.C(uni);
				grid.agregate();

				grid.wyniki();
		
			}

		
	}
	
}

