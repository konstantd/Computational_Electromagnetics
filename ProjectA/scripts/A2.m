%Data
w=3*10^(-2);
t=2*10^(-3);
d=1*10^(-2);
V=100;
er=2.2; 
A=5*w;
e0 = 8.854* 10^(-12); 

%Δημιουργία γεωμετρίας

%Χωρίο Υπολογισμού
rect1 = [3
    4
    -A/2
    -A/2
    +A/2
    +A/2
    -A/2
     A/2
     A/2
    -A/2];

%Κάτω πλάκα
rectdown = [3
    4
    -w/2
    -w/2
    +w/2
    +w/2
    -d/2-t
     -d/2
     -d/2
    -d/2-t];

%Πάνω πλάκα
rectup = [3
    4
    -w/2
    -w/2
    +w/2
    +w/2
    d/2
    d/2+t
    d/2+t
    d/2];

%Xωρίο που περικλειει τις 2 πλακες
rect2 = [3
    4
    -w/2-t
    -w/2-t
    +w/2+t
    +w/2+t
    -d/2-2*t
     d/2+2*t
     d/2+2*t
    -d/2-2*t];

gd = [rect1,rectdown,rectup,rect2];
ns = char('rect1','rectdown','rectup','rect2');
ns = ns';
sf = 'rect1-rectdown-rectup';
d1= decsg(gd,sf,ns);

[p,e,t] = initmesh(d1); 
[p,e,t] = refinemesh(d1,p,e,t); 

pdeplot(p,e,t);     
axis equal

Nn = size(p,2); %Υπολογισμός πλήθους κόμβων
Ne = size(t,2); %Υπολογισμός πλήθους στοιχείων

X0=zeros(Nn,1); %Δηλώνει τις τιμές δυναμικού σε κάθε κόμβο
node_id = ones(Nn,1); %Δηλώνει με 0 αν ο κόμβος ειναι γνωστός, με 1 εάν άγνωστος
Nd = size(e,2);  %Number of edges

for id = 1:Nd  %Τρέχει όλες τις ακμές και ψάχνει τους γνωστούς κόμβους
             
   if  ( (e(6,id)==0 && e(7,id)==2) || ( e(6,id)==2 && e(7,id)==0) )
        if  ( p(2,e(1,id)) > 0 || (p(2,e(2,id)) > 0) )           
            X0(e(1,id))=V/2; 
            X0(e(2,id))=V/2;        
            node_id(e(1,id))=0;  %Mark as known
            node_id(e(2,id))=0;  %Mark as known  
            
         elseif ( (p(2,e(1,id)) < 0)  || (p(2,e(2,id)) < 0) ) %&& (p(2,e(2,id)) == -d/2 ) 
              X0(e(1,id))= -V/2; 
              X0(e(2,id))= -V/2;             
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

for ie = 1:Ne                   %Τρέχει ολα τα στοιχεία
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
            if abs(x(1))<=w/2 && abs(x(2))<=w/2 && abs(x(3))<=w/2 && abs(y(1))<=d/2  && abs(y(2))<=d/2  && abs(y(3))<=d/2
            %if rg==2
                 Se(i,j) = er*e0*(b(i)*b(j)+c(i)*c(j))*Ae;      %Se τοπικός
            else
                 Se(i,j) = e0*(b(i)*b(j)+c(i)*c(j))*Ae;         %Se τοπικός
            end
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

%Βάζουμε και τις τιμές δυναμικού των άγνωστων κόμβων στο Χ0
j=1;
for i=1:Nn
    if (index(i)~=0)             
      X0(i)=X(j);
      j=j+1;
    end
end

%Απεικόνιση Πεδίου
[ux,uy] = pdegrad(p,t,X0); 
pdeplot(p,e,t,'FlowData',[-ux;-uy]);
axis equal
hold on
pdegplot(d1)
axis equal

%Απεικόνιση Δυναμικού
pdeplot(p, e, t, 'XYData', X0);
axis equal

%Υπολογισμός ενέργειας
We=0;
for ie = 1:Ne                   %Τρέχει ολα τα στοιχεία
    n(1:3) = t(1:3,ie);         %3 αριθμοί κόμβων του στοιχείου
    rg = t(4,ie);               %Αριθμός region του στοιχείου 
    x(1:3) = p(1,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    y(1:3) = p(2,n(1:3));       %Συνεταγμένες κόμβων x1, x2, x3
    De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]); 
    Ae = abs(De/2);          %Element Area   
    
    %Συμβάλλουν μονάχα τα στοιχεία ανάμεσα των πλακών
   if abs(x(1))<=w/2 && abs(x(2))<=w/2 && abs(x(3))<=w/2 && abs(y(1))<=d/2 && abs(y(2))<=d/2  && abs(y(3))<=d/2
    We=0.5*er*e0*Ae*(ux(ie)^2 + uy(ie)^2) + We;
   end
end

%Υπολογισμός χωρητικότητας
C=2*We/V^2;

%Ακριβής τιμή χωρητικότητας ανά μονάδα μήκους
Creal=er*e0*w/d;

%Σχετικό Σφάλμα
error=100*(Creal-C)/Creal; 


