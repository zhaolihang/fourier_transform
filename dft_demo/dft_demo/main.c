
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#define PI	3.1415926535


typedef struct				//复数结构体,用于实现傅里叶运算 
{
	double real, image;
} complex;

complex complex_exp_minus(int N, int f, int n)
{
	complex result;
	// 欧拉展开
	double theta = 2 * PI * f * (1 / (double)N) * n;
	result.real = cos(theta);
	result.image = -sin(theta);
	return result;
}

complex complex_exp(int N, int f, int n)
{
	complex result;
	// 欧拉展开
	double theta = 2 * PI * f * (1 / (double)N) * n;
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

typedef struct
{
	int array_count;
	double* array;
	int f_count;
	complex* f_array_result;
	complex* idft_array;
	complex* exp_i;
	complex* exp_minus_i;
} context;

context init_context(int array_count, int f_count) {
	context c;
	c.array_count = array_count;
	c.f_count = f_count;
	c.array = (double*)malloc(sizeof(double) * array_count);
	c.f_array_result = (complex*)malloc(sizeof(complex) * f_count);
	c.idft_array = (complex*)malloc(sizeof(complex) * array_count);

	c.exp_i = (complex*)malloc(sizeof(complex) * array_count);
	c.exp_minus_i = (complex*)malloc(sizeof(complex) * f_count);
	{
		{
			char* flags = malloc(sizeof(char) * array_count);
			memset(flags, 0, sizeof(char) * array_count);
			for (int f = 0; f < f_count; f++)
			{
				for (int n = 0; n < array_count; n++)
				{
					int real_n = f*n % array_count;

					if (flags[real_n] == 0) {
						if (array_count % 2 == 0) {
							int half = array_count / 2;
							if (real_n >= half) {
								int dep_n = real_n - half;
								if (flags[dep_n] == 1) {
									complex result = c.exp_i[dep_n];
									result.real = -result.real;
									result.image = -result.image;
									c.exp_i[real_n] = result;
									flags[n] = 1;
									continue;
								}
							}
						}

						c.exp_i[real_n] = complex_exp(array_count, 1, real_n);
						flags[n] = 1;
					}
				}
			}
			free(flags);
		}
		{
			char* flags = malloc(sizeof(char) * f_count);
			memset(flags, 0, sizeof(char) * f_count);
			for (int f = 0; f < f_count; f++)
			{
				for (int n = 0; n < array_count; n++)
				{
					int real_n = f*n % f_count;

					if (f_count % 2 == 0) {
						int half = f_count / 2;
						if (real_n >= half) {
							int dep_n = real_n - half;
							if (flags[dep_n] == 1) {
								complex result = c.exp_minus_i[dep_n];
								result.real = -result.real;
								result.image = -result.image;
								c.exp_minus_i[real_n] = result;
								flags[n] = 1;
								continue;
							}
						}
					}

					if (flags[real_n] == 0) {
						c.exp_minus_i[real_n] = complex_exp_minus(f_count, 1, real_n);
						flags[n] = 1;
					}
				}
			}
			free(flags);
		}
	}
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
	if (ctx.exp_i != 0) {
		free(ctx.exp_i);
	}
	if (ctx.exp_minus_i != 0) {
		free(ctx.exp_minus_i);
	}
}

complex get_complex_exp_minus(context ctx, int f, int n)
{
	int real_n = f*n % ctx.f_count;
	complex result = ctx.exp_minus_i[real_n];
	return result;
}

complex get_complex_exp(context ctx, int f, int n)
{
	int real_n = f*n % ctx.array_count;
	complex result = ctx.exp_i[real_n];
	return result;
}

complex dft_calculate_f(context ctx,int f)
{
	complex sum;
	sum.real = 0;
	sum.image = 0;
	for (int n = 0; n<ctx.array_count; n++)
	{
		complex result = get_complex_exp_minus(ctx, f, n);
		double origin = ctx.array[n];
		sum.real += origin * result.real;
		sum.image += origin * result.image;
	}
	return sum;
}

void dft_calculate(context ctx)
{
	for (int f = 0; f<ctx.f_count; f++)
	{
		complex result = dft_calculate_f(ctx,f);
		ctx.f_array_result[f] = result;
	}
}

double amplitude(context ctx,int f) {
	complex f_result = ctx.f_array_result[f];
	double real = f_result.real;
	double image = f_result.image;
	double amplitude = sqrt(real * real + image * image);  //计算幅值
	return amplitude;
}

complex idft_calculate_point(context ctx,int n) {
	complex p_result;
	p_result.real = 0;
	p_result.image = 0;

	double _1_N = 1.0 / (double)ctx.array_count;
	for (int f = 0; f<ctx.f_count; f++)
	{
		complex f_result = ctx.f_array_result[f];
		complex complex_ft = get_complex_exp(ctx, f, n);
		complex p_result_tmp = complex_mul(f_result, complex_ft);
		p_result.real += p_result_tmp.real;
		p_result.image += p_result_tmp.image;
	}
	p_result.real *= _1_N;
	p_result.image *= _1_N;
	return p_result;
}

void idft_calculate(context ctx){
	for (int n = 0; n<ctx.array_count; n++)
	{
		complex p_result = idft_calculate_point(ctx,n);
		ctx.idft_array[n] = p_result;
	}
}

int main() {
	int array_count = 16;
	int f_count = array_count;

	context ctx = init_context(array_count, f_count);

	for (int i = 0; i<array_count; i++)//制造输入序列 
	{
		ctx.array[i] = i;
	}

	for (int i = 0; i<array_count; i++)
	{
		printf("%d\t%lf\n", i, ctx.array[i]); //输出计算结果
	}
	printf("===================================================\n");

	dft_calculate(ctx);
	idft_calculate(ctx);

	for (int f = 0; f<f_count; f++)
	{
		complex f_result = ctx.f_array_result[f];
		printf("%d\t%lf\t%lf\t%lf\n", f, f_result.real, f_result.image, amplitude(ctx, f)); //输出计算结果
	}
	printf("===================================================\n");
	for (int i = 0; i<array_count; i++)
	{
		complex p_result = ctx.idft_array[i];
		printf("%d\t%lf\t%lf\n", i, p_result.real, p_result.image); //输出计算结果
	}

	free_context(ctx);
	// matlab 结果对比
	// xn = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]; %输入的采样序列
	// N = 16;  %序列长度
	// f = fft(xn, N);
	// L = 0:1 : N - 1;  %横轴长度
	// subplot(1, 1, 1);  %画图
	// stem(L, abs(f));  %数据源
	// xlabel('k');  %设置横轴名称
	// ylabel('X(k)');   %设置纵轴名称
	// title('DFT N=16');

	getchar();
	return 0;
}