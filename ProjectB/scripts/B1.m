a= 0.01;
newline=sprintf('\n'); 
e0 = 8.854* 10^(-12);  

%���������� ����������
C1 = [1
    0
    0
    a];

gd = [C1];
ns = char('C1');
ns = ns';

sf = 'C1';
dl= decsg(gd,sf,ns);

[p,e,t] = initmesh(dl); 
[p,e,t] = refinemesh(dl,p,e,t); 

pdeplot(p,e,t);     
axis equal

Nn = size(p,2); %����������� ������� ������
Ne = size(t,2); %����������� ������� ���������


X0=zeros(Nn,1); %������� ��� ����� ��������� �� ���� �����
node_id = ones(Nn,1); %������� �� 0 �� � ������ ����� �������, �� 1 ��� ��������
Nd = size(e,2);  %Number of edges

choice = 0;
while ((choice ~= 1) && (choice ~=2))
   
    disp(newline)
    disp('If you want to plot the TE modes press 1.')
    disp('If you want to plot the TM modes press 2.')
    disp(newline)
    prompt = 'Make a choice : ';
    choice = input(prompt);

    while ( choice ~= 1 && choice ~= 2)
        prompt = 'Wrong input! Please try again : ';
        choice = input(prompt);
    end

   if choice == 2
       for id = 1:Nd  %������ ���� ��� ����� ��� ������ ���� �������� �������
          if  e(6,id)==0 || e(7,id)==0      
             node_id(e(1,id))=0;  %Mark as known
             node_id(e(2,id))=0;  %Mark as known        
             X0(e(1,id))=0; 
             X0(e(2,id))=0;   
             end
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


S = spalloc(Nf,Nf,7*Nf); %������ �������
T = spalloc(Nf,Nf,7*Nf); %������ �������

%B = zeros(Nf,1);   %�f �� ������ ��������, 7*Nf ���������� ������ �� ���������

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
            Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;         %Se �������
            if i==j
                Te(i,j)=Ae/6;
            else
                Te(i,j)=Ae/12;
            end
            if (node_id(n(i))~=0)         %� ������ n(i) ����� �������� 
                if (node_id(n(j))~=0)     %� ����� n(j) ����� �������� 
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
                    T(index(n(i)),index(n(j))) = T(index(n(i)),index(n(j))) + Te(i,j);         
                end
            end
        end
    end
end


[V,D] = eigs(S,T,10,'sm');

%����������� 6 ������ ������������ ���������

store=zeros(1,6);   %�� ������� ��� ����� ���� �� ����������������� ��������������� ��� �� ������� ���������
indexx=1;           %�� ������ ��� ���������
counter=0;          %������� ��� ��� ���������� �� ����������������� ��������� 

while indexx<=9
   if  abs(D(indexx,indexx)-D(indexx+1,indexx+1) )/10^5 <1  
      counter = counter+1;                   %���������� ��� �������
      kc(counter)=sqrt(D(indexx,indexx));   %������� ��������� ������ �������� 
      store(counter)=indexx;    %���������� ������
      indexx=indexx+2;          %������ ��� ������� ��������
  else 
     counter = counter+1;
     kc(counter)=sqrt(D(indexx,indexx));
     store(counter)=indexx;
     indexx=indexx+1;          
end
end
 
%����� � ��������� ��������
kc(6)=sqrt(D(10,10));
store(6)=10;

%���������� ��������
fc=zeros(1,6);    
fc=3*10^8*kc/(2*pi);

%����������� ��� �������� ���������
if choice==2
   pnm_TM=[7.016 6.38 5.52 5.135 3.832 2.405];
   fcreal_TM=3*(10^8)*pnm_TM/(2*pi*a);      %���������� ���������� ��������
    error_TM=100*abs((fcreal_TM-fc))./fcreal_TM;
else
    pnm_TE=[5.317 4.201 3.832 3.054 1.841 0];
    fcreal_TE=3*(10^8)*pnm_TE/(2*pi*a);
    error_TE=100*abs((fcreal_TE-fc))./fcreal_TE;

end



for q=1:6 %��� 6 �������

%������� ��� ��� ����� ��������� ��� �������� ������ ��� �0
j=1;
for i=1:Nn
    if (index(i)~=0) %if ��������      
      X0(i)=V(j,store(q));
      j=j+1;
    end
end

%���������� ������
figure
pdeplot(p, e, t, 'XYData', X0);
axis equal tight
colormap jet

end




