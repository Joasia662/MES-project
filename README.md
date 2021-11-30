# MES-project

The Main goal of the program is to determine the temperature distribution in the cross-section - maximum and minimum temperature at
set initial temperature, ambient temperature and parameters of a given body (including specific heat, density, etc.).

<b> Description of the created classes and their role: </b>
* Node: reflects a single point, contains x, y, t, node id and flag which tells us whether the knot is near the shore or not.
* UniversalEl: Used to compute and store shape functions and their derivatives above ksi and eta in the universal element.
* Point: integration points, contains values and weights.
* Grid: grid of elements, it is used to create all elements, and then to aggregation of matrices and vectors to global.
* Globals: static class for storing global data
* Function_of_Shape - static class that stores shape functions
* Main: we start our program from here, we call the mesh which further creates elements and calls calculation functions
* GaussianElimination: class containing a mathematical algorithm for solving systems linear equations, calculating the matrix momentum, calculating the inverse matrix 


------------------------------------------------------------------------
Celem programu jest wyznaczenie rozkładu temperatury w przekroju- maksymalnej i minimalnej temperatury przy
zadanej temperaturze początkowej, temperaturze otoczenia oraz parametrach danego ciała (w tymciepła właściwego, gęstości itp.).

<b> Opis stworzonych klas oraz ich rola:</b>
* Node: odzwierciedla pojedyńczy punkt, zawiera rzeczywiste x,y,t, id węzła oraz flagę która mówi nam czy węzęł jest przy brzegu czy nie.
* UniversalEl: Służy do wyliczania,przechowywania funkcji kształtu i ich pochodnych po ksi oraz eta w elemencie uniwersalym.
* Point: czyli punkty całkowania, zawiera wartości oraz wagi.
* Siatka: siatka elementów, służy do stworzenia wszystkich elementów, a potem do agregacji macierzy oraz wektorów na globalny.
* Globals : statyczna klasa do przechowywania danych globalnych
* Function_of_Shape - statyczna klasa przechowująca funkcje kształtu
* Main: startujemy stąd nasz program, wywołujemy siatkę która dalej tworzy elementy oraz wywołuje funkcje obliczeń
* GaussianElimination: klasa zawierająca matematyczny algorytm rozwiązywania układów równań liniowych, obliczania zędu macierzy, obliczania macierzy odwrotnej

Metoda Elementów Skończonych jest metodą aproksymacji (czyli otrzymywania rozwiązań przybliżonych) równań różniczkowych cząstkowych (RRC). Równania różniczkowe stanowią model matematyczny,najczęściej jakiegoś procesu lub stanu układu fizycznego. Proces lub stan opisywane są za pomocą parametrów będących funkcjami położenia w przestrzeni i ewentualnie czasu. <b> Zastosowanie MES do rozwiązania konkretnego zadania naukowego lub inżynierskiego składa się z dwóch odrębnych procesów:
stworzenia modelu obliczeniowego rozwiązania konkretnego zadania za pomocą uzyskanego modelu. </b>
