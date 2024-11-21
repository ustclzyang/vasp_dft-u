% 计算U值
t=readmatrix("mylog");
t=sort(t);
x=t(:,1);
y1=t(:,2);
y2=t(:,3);
p1=polyfit(x,y1,1);
p2=polyfit(x,y2,1);
figure
plot(x,y1,'DisplayName',sprintf('NSCF: Y=%.6f N_d + %.3f',p1(1),p1(2)),'Color','k','Marker','o')
hold on
plot(x,y2,'DisplayName',sprintf('SCF: Y=%.6f N_d + %.3f',p2(1),p2(2)),'Color','r','Marker','square')
legend('Location','northwest')
ylabel('Number of d-electrons')
xlabel('V (eV)')
U=1/p2(1)-1/p1(1);
