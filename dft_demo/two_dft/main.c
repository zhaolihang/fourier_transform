
#include <stdio.h>
#include <stdlib.h>
#include <math.h>			 
#define PI	3.1415926535


typedef struct				//复数结构体,用于实现傅里叶运算 
{
	double real, image;
} complex;


typedef struct
{
	int array_x_count;
	int array_y_count;

	double* array;
	int f_x_count;
	int f_y_count;
	complex* f_array_result;
	complex* idft_array;
} context;

context init_context(int array_x_count, int array_y_count, int f_x_count, int f_y_count) {
	context c;
	c.array_x_count = array_x_count;
	c.array_y_count = array_y_count;
	c.f_x_count = f_x_count;
	c.f_y_count = f_y_count;
	c.array = (double*)malloc(sizeof(double) * array_x_count*array_y_count);
	c.f_array_result = (complex*)malloc(sizeof(complex) * f_x_count*f_y_count);
	c.idft_array = (complex*)malloc(sizeof(complex) * array_x_count*array_y_count);
	return c;
}

void free_context(context ctx) {
	if (ctx.array != 0) {
		free(ctx.array);
	}
	if (ctx.f_array_result != 0) {
		free(ctx.f_array_result);
	}
	if (ctx.idft_array != 0) {
		free(ctx.idft_array);
	}
}

complex complex_exp_minus(int fx,int Nx ,int nx, int fy, int Ny, int ny)
{
	complex result;
	// 欧拉展开
	double theta = 2 * PI * (fx * (1 / (double)Nx) * nx + fy * (1 / (double)Ny) * ny);
	result.real = cos(theta);
	result.image = -sin(theta);
	return result;
}

complex complex_exp(int fx, int Nx, int nx, int fy, int Ny, int ny)
{
	complex result;
	// 欧拉展开
	double theta = 2 * PI * (fx * (1 / (double)Nx) * nx + fy * (1 / (double)Ny) * ny);
	result.real = cos(theta);
	result.image = sin(theta);
	return result;
}

complex complex_mul(complex a, complex b)
{
	complex result;
	result.real = a.real * b.real - a.image * b.image;
	result.image = a.real * b.image + a.image * b.real;
	return result;
}

int get_array_index(context ctx,int nx,int ny)
{
	return ctx.array_y_count * nx + ny;
}

int get_f_array_result_index(context ctx, int fx, int fy)
{
	return ctx.f_y_count * fx + fy;
}

complex dft_calculate_f(context ctx,int fx,int fy)
{
	complex sum;
	sum.real = 0;
	sum.image = 0;

	for (int nx = 0; nx<ctx.array_x_count; nx++)
	{
		for (int ny = 0; ny<ctx.array_y_count; ny++)
		{
			complex result = complex_exp_minus(fx, ctx.array_x_count, nx, fy, ctx.array_y_count, ny);
			double origin = ctx.array[get_array_index(ctx, nx, ny)];
			sum.real += origin * result.real;
			sum.image += origin * result.image;
		}
	}
	return sum;
}

void dft_calculate(context ctx)
{
	for (int fx = 0; fx<ctx.f_x_count; fx++)
	{
		for (int fy = 0; fy<ctx.f_y_count; fy++)
		{
			complex result = dft_calculate_f(ctx, fx, fy);
			ctx.f_array_result[get_f_array_result_index(ctx, fx, fy)] = result;
		}
	}
}

double amplitude(context ctx, int fx, int fy) {
	complex f_result = ctx.f_array_result[get_f_array_result_index(ctx, fx, fy)];
	double real = f_result.real;
	double image = f_result.image;
	double amplitude = sqrt(real * real + image * image);  //计算幅值
	return amplitude;
}


complex idft_calculate_point(context ctx, int nx, int ny)
{
	complex p_result;
	p_result.real = 0;
	p_result.image = 0;

	double _1_Nx = 1.0 / (double)ctx.array_x_count;
	double _1_Ny = 1.0 / (double)ctx.array_y_count;

	for (int fx = 0; fx<ctx.f_x_count; fx++)
	{
		for (int fy = 0; fy<ctx.f_y_count; fy++)
		{
			complex f_result = ctx.f_array_result[get_f_array_result_index(ctx,fx,fy)];
			complex complex_ft = complex_exp(fx, ctx.array_x_count, nx, fy, ctx.array_y_count, ny);
			complex p_result_tmp = complex_mul(f_result, complex_ft);
			p_result.real += p_result_tmp.real;
			p_result.image += p_result_tmp.image;
		}
	}

	p_result.real *= _1_Nx;
	p_result.image *= _1_Nx;
	p_result.real *= _1_Ny;
	p_result.image *= _1_Ny;
	return p_result;
}

void idft_calculate(context ctx) {
	for (int nx = 0; nx<ctx.array_x_count; nx++)
	{
		for (int ny = 0; ny<ctx.array_y_count; ny++)
		{
			complex p_result = idft_calculate_point(ctx, nx, ny);
			ctx.idft_array[get_array_index(ctx, nx, ny)] = p_result;
		}
	}
}

int main() {

	int array_x_count = 4;
	int array_y_count = 4;

	int f_x_count = array_x_count;
	int f_y_count = array_y_count;

	context ctx = init_context(array_x_count, array_y_count, f_x_count, f_y_count);
	for (int nx = 0; nx<ctx.array_x_count; nx++)
	{
		for (int ny = 0; ny<ctx.array_y_count; ny++)
		{
			ctx.array[get_array_index(ctx, nx, ny)] = get_array_index(ctx, nx, ny);
			printf("%d=%lf\n", get_array_index(ctx, nx, ny), ctx.array[get_array_index(ctx, nx, ny)]);
		}
	}
	printf("===================================================\n");

	dft_calculate(ctx);
	idft_calculate(ctx);

	for (int fx = 0; fx<ctx.f_x_count; fx++)
	{
		for (int fy = 0; fy<ctx.f_y_count; fy++)
		{
			complex f_result = ctx.f_array_result[get_f_array_result_index(ctx, fx, fy)];
			printf("%d,%d\t%lf\t%lf\t%lf\n", fx, fy, f_result.real, f_result.image, amplitude(ctx, fx, fy)); //输出计算结果
		}
	}

	printf("===================================================\n");

	for (int nx = 0; nx<ctx.array_x_count; nx++)
	{
		for (int ny = 0; ny<ctx.array_y_count; ny++)
		{
			int array_index = get_array_index(ctx, nx, ny);
			complex p_result = ctx.idft_array[array_index];
			printf("%d,%d\t%lf\t%lf\n", nx, ny, p_result.real, p_result.image); //输出计算结果
		}
	}

	free_context(ctx);
	// matlab 结果对比
	// xn = [0, 1, 2, 3; 4, 5, 6, 7; 8, 9, 10, 11; 12, 13, 14, 15];
	// Nx = 4; Ny = 4;
	// f = fft2(xn, Nx, Ny);
	// figure(1);
	// imshow(xn, [0 15], 'InitialMagnification', 'fit'); %显示图像
	// figure(2);
	// F2 = log(abs(f)); %对傅里叶变换结果取对数可视化
	// imshow(F2, [0 5], 'InitialMagnification', 'fit'); %显示图像

	getchar();
	return 0;
}