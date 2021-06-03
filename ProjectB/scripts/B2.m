
 %Data
lambda=1;
a_values=[lambda/4 lambda 5*lambda/2];
d_values=[lambda lambda/4 ];
 
 for a_index=1:size(a_values,2)
   for d_index=1:size(d_values,2)

f=300e+6;
d=d_values(d_index);
a=a_values(a_index);
Eo=1;
w=1;
b1=log(10^-6)/(-2*d*2*pi);
a1=complex(1,-b1);

e0 = 8.854187817* 10^(-12);  
m0=1.256637061*10^(-6);

%���������� ����������

%�������� ���������
C1 = [1
    0
    0
    a];

%����������� ������� ��
rect2 = [3
    4
    -a-w
    -a-w
    a+w
    a+w
    -a-w
    a+w
    a+w
    -a-w];


%�������� PML

%�������� �x
rect3 = [3
    4
    -a-w-d
    -a-w-d
    -a-w
    -a-w
    -a-w
    a+w
    a+w
    -a-w];

%���� �������� �x�y
rect4 = [3
    4
    -a-w-d
    -a-w-d
    -a-w
    -a-w
    a+w
    a+w+d
    a+w+d
    a+w];

%���� �������� �x�y
rect5 = [3
    4
    -a-w-d
    -a-w-d
    -a-w
    -a-w
     -a-w-d
     -a-w
      -a-w
   -a-w-d];

%����� �x
rect6 = [3
    4
 a+w
 a+w
 a+w+d
 a+w+d
 -a-w
 a+w
 a+w
 -a-w];


%���� ����� �x�y
rect7 = [3
    4
 a+w
 a+w
 a+w+d
 a+w+d
a+w
a+w+d
+a+w+d
a+w];

%���� ����� �x�y
rect8 = [3
    4
 a+w
 a+w
 a+w+d
 a+w+d
-a-w-d
-a-w
-a-w
-a-w-d];

%���� �y
rect9=[
    3
    4
    -a-w
    -a-w
    a+w
    a+w
    a+w
    a+w+d
    a+w+d
    a+w];

%���� �y
rect10=[
    3
    4
    -a-w
    -a-w
    a+w
    a+w
   -a-d-w
   -a-w
   -a-w
   -a-w-d];

C1 = [C1;zeros(length(rect2) - length(C1),1)];  %���� ������� �� �� �������� �������
gd = [C1, rect2, rect3, rect4, rect5, rect6, rect7, rect8, rect9, rect10];
ns = char('C1','rect2','rect3','rect4','rect5','rect6','rect7','rect8', 'rect9', 'rect10');
ns = ns';

sf = '(rect2+rect3+rect4+rect5+rect6+rect7+rect8+rect9+rect10)-C1';
dl= decsg(gd,sf,ns);

[p,e,t] = initmesh(dl); 
[p,e,t] = refinemesh(dl,p,e,t); 
%[p,e,t] = refinemesh(dl,p,e,t); 
%[p,e,t] = refinemesh(dl,p,e,t); 
%[p,e,t] = refinemesh(dl,p,e,t); 

Nn = size(p,2); %����������� ������� ������
Ne = size(t,2); %����������� ������� ���������


X0=zeros(Nn,1); %������� ��� ����� ��������� �� ���� �����
node_id = ones(Nn,1); %������� �� 0 �� � ������ ����� �������, �� 1 ��� ��������
Nd = size(e,2);  %Number of edges



for id = 1:Nd  %������ ���� ��� ����� ��� ������ ���� �������� �������
    if  e(6,id)==0 || e(7,id)==0      
       if   sqrt( p(1,e(1,id))^2 + p(2,e(1,id))^2 )==a
       
        node_id(e(1,id))=0;  %Mark as known
        node_id(e(2,id))=0;  %Mark as known
        
        X0(e(1,id))= -Eo*exp(-1i*2*pi*p(1,e(1,id))); 
        X0(e(2,id))= -Eo*exp(-1i*2*pi*p(1,e(2,id)));         
       
    end 
    end
end


index=zeros(Nn,1); %������� ��� ������ ��� �������� ������ ��� node_id/X0
                   %������� ������ �� ����� index �����
counter = 0; 
for in= 1:Nn   %������ ����� ���� �������    
      if node_id(in) == 1   %��������
         counter = counter + 1;
        index(in)=counter; 
    end
end

Nf=counter; %������ ��������
A = spalloc(Nf,Nf,7*Nf); %������ ������� 
B= zeros(Nf,1);   %�f �� ������ ��������, 7*Nf ���������� ������ �� ���������


for ie = 1:Ne                  %������ ��� �� ��������
    n(1:3) = t(1:3,ie);         %3 ������� ������ ��� ���������
    rg = t(4,ie);               %������� region ��� ��������� 
    x(1:3) = p(1,n(1:3));       %������������ ������ x1, x2, x3
    y(1:3) = p(2,n(1:3));       %������������ ������ x1, x2, x3
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
           
             if rg==5  %�����
                    ezz=e0;
                    mxx=m0;
                    myy=m0;
                  
                elseif rg==3 || rg==8 || rg==1 || rg==9  %�x�y
                    ezz= a1^2*e0;
                    mxx=m0;
                    myy=m0;
                   
                elseif rg==4 ||  rg==6 %�x 
                     ezz= a1*e0;
                     mxx=m0/a1;
                     myy=m0*a1;
           
                elseif rg==2 || rg==7 %�y
                      ezz= a1*e0;
                      myy=m0/a1;            
                      mxx=a1*m0;
               end
        Se(i,j) = (myy^(-1)*b(i)*b(j)+mxx^(-1)*c(i)*c(j))*Ae; 
        
             if i==j
                Te(i,j)=ezz*Ae/6;
            else
                Te(i,j)=ezz*Ae/12;
            end
         AE(i,j)= Se(i,j)-((2*pi*f)^2)*Te(i,j);
           
            if (node_id(n(i))~=0)         %� ������ n(i) ����� �������� 
                if (node_id(n(j))~=0)     %� ����� n(j) ����� �������� 
                    A(index(n(i)),index(n(j))) = A(index(n(i)),index(n(j))) + AE(i,j) ;
                else
                    B(index(n(i))) = B(index(n(i)))-AE(i,j)*X0(n(j)); 
                end
            end
        end
    end
end

%���� ����������
Es=A\B;

%����������� �� ���� ����� �� ���������� 
for i=1:Nn
 
    if (index(i)~=0)   %if ��������
      X0(i)= abs(exp(-1i*2*pi*p(1,i))+ Es(index(i)));
    else
      X0(i)= abs( X0(i) + exp(-1i*2*pi*p(1,i)));  
    end

end


figure
pdeplot(p, e, t, 'XYData', X0);
axis tight
axis equal
%��� ����������������� ��� PML-��������
xlim([-a-w   a+w])
ylim([-a-w   a+w])
str = sprintf('d = %1.2f ,a = %1.2f', d, a);
title(str)
colormap jet



   end
end

