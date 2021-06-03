%Data
b=3.5*10^(-3) / 2;
Zo=50;
a=b/(exp(Zo/60));
V=1;
e0 = 8.854* 10^(-12);  

%Δημιουργία γεωμετρίας

%Εσωτερικός αγωγός
C1 = [1
    0
    0
    a];

%Εξωτερικός αγωγός
C2 = [1
    0
    0
    b];

gd = [C1,C2];
ns = char('C1','C2');
ns = ns';

sf = 'C2-C1';
dl= decsg(gd,sf,ns);

[p,e,t] = initmesh(dl); 
[p,e,t] = refinemesh(dl,p,e,t); 
[p,e,t] = refinemesh(dl,p,e,t); 
[p,e,t] = refinemesh(dl,p,e,t); 

pdeplot(p,e,t);     
axis equal

Nn = size(p,2); % Υπολογισμός πλήθους κόμβων
Ne = size(t,2); %Υπολογισμός πλήθους στοιχείων

X0=zeros(Nn,1); %Δηλώνει τις τιμές δυναμικού σε κάθε κόμβο
node_id = ones(Nn,1); %Δηλώνει με 0 αν ο κόμβος ειναι γνωστός, με 1 εάν άγνωστος
Nd = size(e,2);  %Number of edges

for id = 1:Nd  %Τρέχει όλες τις ακμές και ψάχνει τους γνωστούς κόμβους
    if  e(6,id)==0 || e(7,id)==0      
       if   sqrt( p(1,e(1,id))^2 + p(2,e(1,id))^2 )==a
        node_id(e(1,id))=0;  %Mark as known
        node_id(e(2,id))=0;  %Mark as known
        X0(e(1,id))=V; 
        X0(e(2,id))=V;         
       elseif  sqrt( p(1,e(1,id))^2 + p(2,e(1,id))^2 )==b  %Ολοι οι υπολοιποι προς το παρον 0
        node_id(e(1,id))=0;  %Mark as known
        node_id(e(2,id))=0;  %Mark as known
      end
    end 
end

index=zeros(Nn,1); %Κρατάει τις θέσεις των άγνωστων κόμβων στο node_id/X0
                  %Γνωστοί κόμβοι θα έχουν index μηδέν
counter = 0; 
for in= 1:Nn   %Τρέχει όλους τους κόμβους    
      if node_id(in) == 1   %’γνωστος
         counter = counter + 1;
        index(in)=counter; 
    end
end

Nf=counter; %Πλήθος αγνώστων
S = spalloc(Nf,Nf,7*Nf); %Αραιός πίνακας
B = zeros(Nf,1);   %Νf το πλήθος αγνώστων, 7*Nf εκτιμώμενο πλήθος μη μηδενικών

for ie = 1:Ne                  %Τρέχει ολα τα στοιχεία
    n(1:3) = t(1:3,ie);         %3 αριθμοί κόμβων του στοιχείου
    rg = t(4,ie);               %Αριθμός region του στοιχείου 
    x(1:3) = p(1,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    y(1:3) = p(2,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]); 
    Ae = abs(De/2);          %Element Area
    b(1) = (y(2)-y(3))/De;
    c(1) = (x(3)-x(2))/De;
    b(2) = (y(3)-y(1))/De;
    c(2) = (x(1)-x(3))/De; 
    b(3) = (y(1)-y(2))/De;     
    c(3) = (x(2)-x(1))/De;   
    for i=1:3 
        for j=1:3 
            Se(i,j) = e0*(b(i)*b(j)+c(i)*c(j))*Ae;         %Se τοπικός
            if (node_id(n(i))~=0)         %Ο κόμβος n(i) είναι άγνωστος 
                if (node_id(n(j))~=0)     %Ο κόβος n(j) είναι άγνωστος 
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
                else
                    B(index(n(i))) = B(index(n(i)))-Se(i,j)*X0(n(j)); 
                end
            end
        end
    end
end


%Λύση Συστήματος
X=S\B;

tic
%X = bicg(S,B,1e-6,100);
toc

%Βάζουμε και τις τιμές δυναμικού των άγνωστων κόμβων στο Χ0
j=1;
for i=1:Nn
    if (index(i)~=0) %if ’γνωστος      
      X0(i)=X(j);
      j=j+1;
    end
end

%Απεικόνιση Πεδίου
[ux,uy] = pdegrad(p,t,X0); 
pdeplot(p,e,t,'FlowData',[-ux;-uy]);
axis equal
hold on
pdegplot(dl)
hold off


%Απεικόνιση Δυναμικού
pdeplot(p, e, t, 'XYData', X0);
axis equal


%Ξαναδηλώνουμε το b, καθώς χρησιμοποιήθηκε στην επίλυση του συστήματος
b=3.5*10^(-3) / 2;

%Υπολογισμός ενέργειας
We=0;
for ie = 1:Ne                  %Τρέχει ολα τα στοιχεία
    n(1:3) = t(1:3,ie);         %3 αριθμοί κόμβων του στοιχείου
    rg = t(4,ie);               %Αριθμός region του στοιχείου 
    x(1:3) = p(1,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    y(1:3) = p(2,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]); 
    Ae = abs(De/2);          %Element Area   
    
   
    We=e0*Ae*(ux(ie)^2 + uy(ie)^2) +We;
end

%Υπολογισμός χωρητικότητας
C=We; %Αφου δεν πολλαπλασίασα με 0.5

%Ακριβής τιμή χωρητικότητας ανά μονάδα μήκους
Creal=2*pi*e0/(log(b/a));

%Σχετικό Σφάλμα
error=100*(Creal-C)/Creal; 


%pdeplot(p, e, t, 'XYData', X0);
axis equal
