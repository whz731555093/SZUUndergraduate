%% ��λ����1��ȡ�������ֱ������ƽ��
clear all;clc;
bx1=0; by1=0;  
bx1=-5e-9; by1=-5e-9;    
bx2=5e-9; by2=5e-9;     
cx=3e-10; cy=3e-10;
d=4e-9;
N=250;
dx0=0.096e-9; dy0=0.096e-9;

xx1 = (-N/2:N/2-1) * dx0; % Input f0 is in natural order
yy1 = (-N/2:N/2-1) * dy0;
[x0,y0]=meshgrid(xx1,yy1);
%y1=exp(-((x0-bx1).^2/(2*cx^2)+(y0-by1).^2/(2*cy^2))) + exp(-((x0-bx2).^2/(2*cx^2)+(y0-by2).^2/(2*cy^2)));
y1=1/(2*pi*cx*cy)*exp(-((x0-bx1).^2/(2*cx^2)+(y0-by1).^2/(2*cy^2))) + 1/(2*pi*cx*cy)*exp(-((x0-bx2).^2/(2*cx^2)+(y0-by2).^2/(2*cy^2)));
figure(1);imshow(y1);

y2=zeros(N,N);
y2(:,N/2)=exp(-xx1/d);
f=conv2(y1,y2);
tj=1;dd=zeros(4,tj);aaa=1;bbb=1;ccc=aaa+bbb+1;dnn=1;
%imshow(mat2gray(f))

A1 = imnoise(mat2gray(f),'poisson');
A2 = imnoise(mat2gray(A1),'gaussian');
A3=wiener2(A2,[3 3]); %step1,filter
%A3 = SF3(1:N,N/2:N+N/2-1);
imshow(A3);
%% cftool���x���˹����
imagex = sum(A3,1);   %��ͼ���������ˮƽ����x�ᣩͶӰ����
res = double(imagex);   % ������ת��Ϊdouble����
[cols, rows] = size(res);   % �������ݾ���ߴ磬�˴����ݼ�Ϊ1��120��
gsx_axis = 1:size(res,2);   % X�����������
gsy_axis = res/max(res);     % Y��Ϊ���ݵ��ֵ
cftool
%% cftool���y���˹����
imagey = sum(y1,2);   %��ͼ���������ˮƽ����x�ᣩͶӰ����
res = double(imagey);   % ������ת��Ϊdouble����
[cols, rows] = size(res);   % �������ݾ���ߴ磬�˴����ݼ�Ϊ1��120��
gsx_axis = 1:size(res,1);   % X�����������
gsy_axis = res/max(res);     % Y��Ϊ���ݵ��ֵ
cftool
%% lsqcurvefit��ϸ�˹�����������˥��������ʱ�䳣���ĳ�ֵ
y1=73.92; y2=178.1;
sigmay1=4.419; sigmay2=4.419;
x1=198.0; x2=302.0;
sigmax1=4.118; sigmax2=4.237;
a=round(x1);
b=round(x2);
imagee=A3(:,a)./6+A3(:,a-1)./6+A3(:,a+1)./6+A3(:,b)./6+A3(:,b-1)./6+A3(:,b+1)./6;   %ȡ������н�����϶�ʱ�䳣����ֵ
res = double(imagee);   % ������ת��Ϊdouble����
res_y = res';   %��������ת��Ϊ������
[rows, cols] = size(res_y);   % �������ݾ���ߴ磬�˴����ݼ�Ϊ1��239��  rows=1��cols=239
ex_axis = 1:size(res_y,2);   %x����������
ey_axis = res_y/max(res_y);   %y�����ݹ�һ��

oy(1)=y1;   %��˹������ֵΪy1=[52.2527103069531]
oy(2)=y2;   %��˹������ֵΪy2=[61.0450913412387]
alocation = find(res_y==max(res_y));   %�������ֵ��������ֵͬʱ���ڵ����������λ��
halfpeaks=find(abs(res_y(:,1:alocation(1))-max(res_y)./2)==min(abs(res_y(:,1:alocation(1))-max(res_y)./2)));  %�ҵ���ߵĴ���˹�İ��ֵλ��=57
elocation=find(abs(res_y(:,alocation(1):cols)-max(res_y)./2.718281)==min(abs(res_y(:,alocation(1):cols)-max(res_y)./2.718281)));  %%����ʱ�䳣�������ͬʱ���ڵ����������λ��
oy(3)=elocation(1);  %ӫ������������ֵ=[13.4586529387669]
oy(4)=abs(alocation(1)-halfpeaks(1))*2/2.35482;   %��׼��=[8.49321816529501];cftool��ϵı�׼��sigmay1=2.249
oy(5)=0;   %����=[0.156884137631035]
oy(6)=1;   %y1�����ֵ=[5.85843486636557]
oy(7)=1;   %y2�����ֵ=[-3.41001445017942]
fun2=@(oy,x)oy(5)+oy(6)*conv(exp(-((x-oy(1)).^2)./(2*oy(3).^2)),exp(-x./oy(3)),'same')+...
    oy(7)*conv(exp(-((x-oy(2)).^2)./(2*oy(3).^2)),exp(-x./oy(3)),'same');%�Ը�˹����*��ָ��˥����������ĺ��������
options=optimset('TolX',2.000000e+02); 
[yy,resnormy] = lsqcurvefit(fun2,oy,ex_axis,ey_axis,[],[],options);
yy
%% �������
facex_axis = 1:size(A3,2);   %ȡx���������䣨����120��
facey_axis = 1:size(A3,1);   %ȡy���������䣨����239��
[X,Y] = meshgrid(facex_axis,facey_axis);
facez_axis = double(A3);    %��image3����ת��Ϊdouble����
xy=[X(:) Y(:)];
z=facez_axis(:);

oi(1)=x1;   %x1=56.8499959303641
oi(2)=x2;   %x2=56.8500000095986
oi(3)=yy(1);   %y1=51.9810100467012
oi(4)=yy(2);   %y2=60.7482279584067
oi(5)=yy(3);   %ӫ������=13.4586529361055
oi(6)=(sigmax1+sigmay1)/2;   %x1y1��sigma1=2.64248398340200
oi(7)=(sigmax2+sigmay2)/2;  %x2y2��sigma2=3.04349995829075
oi(8)=1;   %���=1.02744649527525
oi(9)=1;   %���=0.954276775324733
oi(10)=0;   %����=0
fun3=@(oi,x) oi(10)+conv(oi(8).*exp(-((x(:,1)-oi(1)).^2+(x(:,2)-oi(3).^2))./(2.*oi(6).^2)),exp(-x(:,2)./oi(5)),'same')+...
    conv(oi(9).*exp(-((x(:,1)-oi(2)).^2+(x(:,2)-oi(4).^2))./(2.*oi(7).^2)),exp(-x(:,2)./oi(5)),'same');
options=optimset('TolX',2.014195e-01); 
[ii,resonrm]=lsqcurvefit(fun3,oi,xy,z,[],[],options);   %resonrm=[5.47346865789654e+179]
ii
%% fitting
sf1=sum(A3);  %A3��image3��sf1��imagex
figure(2);subplot(1,2,1);imshow(A2);subplot(1,2,2);imshow(SF3);
rm=find(sf1==max(sf1));   %��image3�ķ�ֵ

cm0(1)=xx1(rm);   %����������ֵ�ĺ����괢����cm0��1��
rrm=find(abs(sf1(1:N/2)-max(sf1)/2)==min(abs(sf1(1:N/2)-max(sf1)/2)));   %�ҵ����ֵ��λ��
cm0(2)=abs(cm0(1)-xx1(rrm))*2/2.35482;%%%��׼��
cm0(3)=0;

fun1 = @(cm,x) cm(3)+exp(-(x-cm(1)).^2/(2*cm(2)^2));
options=optimset('TolX',1e-10^5); 
ccm=lsqcurvefit(fun1,[cm0(1) cm0(2) cm0(3)],xx1,sf1,[],[],options);

rmm=find(abs(xx1(1:N/2)-ccm(1))==min(abs(xx1(1:N/2)-ccm(1))));%%rrmɶ��˼

%����˥�����
t=1;   cc=zeros(5,ccc);   ee=0;
for i=N/2+rmm-aaa:N/2+rmm+bbb  %%%%%��һ��ɶ��˼
Aa=SF3(:,i);   %��i�������й��ɵ�������
A=Aa';  %����i������תΪ������
% figure(3);subplot(1,ccc,t);plot(Aa,'k');hold on;
B=A(1:N);   %��������A�ĵ�1����N����Ԫ����ɵ���������

c0(1)=1;   %start-guess here 
r1=find(B==max(B));   %�ҵ�y�����ֵλ�ü�¼��r1,B=res_y
c0(2)=xx1(r1(1));  %����ֵ��ֵ������c0(2),y(1)Ϊ��ֵ��ֵ
rr=find(abs(B(1:r1(1))-max(B)/2)==min(abs(B(1:r1(1))-max(B)/2)));  %����߰�ߣ��������ҵ���ߵĴ���˹�İ��ֵ��
rrr=find(abs(B(r1(1):N)-max(B)/2.71828)==min(abs(B(r1(1):N)-max(B)/2.71828)));  %�ұ߰��
c0(3)=abs(c0(2)-xx1(rr))*2/2.35482;  %��׼��
c0(4)=abs(xx1(rrr+r1(1)-1)-c0(2));  %����
c0(5)=0;   %����

fun = @(c,x) c(5)+conv(c(1)*exp(-(x-c(2)).^2/(2*c(3)^2)),exp(-x/c(4)));   %�Ը�˹����*��ָ��˥����������ĺ��������
options=optimset('TolX',1e-10^5); 
[cc,res]=lsqcurvefit(fun,c0,xx1,A,[],[],options);

if (res<1)
    ee=cc+ee;
    t=t+1;%%%���ѭ��ɶ��˼
end
    Ifit=fun(cc,xx1); %your fitted gaussian in vector 
%     plot(Ifit,'r')
end
% ddcc(dn,1)=res;ddcc(dn,2:4)=cc(2:4);
ddj=ee/(t-1);
if (t>1)
    dd(1:3,dnn)=ddj(2:4);
    dd(4,dnn)=ccm(1);
    dnn=dnn+1;%%%%���ѭ��ɶ��˼
end

%end


ddd=dd';
xlswrite('dd.xls',ddd);
ds=find(dd(1,:)==0);
if numel(ds)==0%%%����Ԫ�ظ���
    ds(1)=tj+1;
end
dds(1:(ds(1)-1),1:4)=ddd(1:(ds(1)-1),1:4);
figure()
[locx1,locx2]=hist(dds(:,4),20);%localization accuracy in x
dealfitx(locx2,locx1);
figure()
[locy1,locy2]=hist(dds(:,1),20);%localization accuracy in y
dealfity(locy2,locy1);
% [NN1,XX1]=hist(dd(1,:),10);
figure()
[loct1,loct2]=hist(dds(:,3),20);title('tao= 1000e-9')%localization accuracy in tao
dealfitt(loct2,loct1);
% [NN3,XX3]=hist(dd(3,:),10)
cftool