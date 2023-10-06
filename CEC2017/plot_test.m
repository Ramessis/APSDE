[x,y] = meshgrid(-5:0.5:5); % 快速生成网格所需的数据
tem = sqrt(x.^2+y.^2)+1e-12; 
z = sin(tem)./tem; % 如果不对tem处理，那么z的最中间的一个值 0/0 = NaN
subplot(1,2,1)
mesh(x,y,z)
xlabel('x轴'); ylabel('y轴'); zlabel('z轴'); % 加上坐标轴的标签
axis vis3d % 冻结屏幕高宽比，使得一个三维对象的旋转不会改变坐标轴的刻度显示
title('mesh(x,y,z)')
subplot(1,2,2)
surf(x,y,z) % (X(j), Y(i), Z(i,j))是线框网格线的交点
xlabel('x轴'); ylabel('y轴'); zlabel('z轴'); % 加上坐标轴的标签
axis vis3d % 冻结屏幕高宽比，使得一个三维对象的旋转不会改变坐标轴的刻度显示
title('surf(x,y,z)')