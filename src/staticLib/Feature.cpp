#include "defines.h"
#include "Feature.h"

Feature::Feature(void) {
	x = 0;
	y = 0;
	scale = 0;
	orientation = 0;
	peak = 0;
	curvRatio = 0;
	pdf = 0;
	spatialIndex = 0;
	relevance = 0;
	octave = 0;
	iscale = 0;
	descr[0] = -1.0f;		// trap to detect uninitialized descriptors when compressing
}

void Feature::toFile(FILE * file) const {
	size_t fout;
	fout = fwrite(&x, sizeof(x), 1, file);
	fout = fwrite(&y, sizeof(y), 1, file);
	fout = fwrite(&scale, sizeof(scale), 1, file);
	fout = fwrite(&orientation, sizeof(orientation), 1, file);
	fout = fwrite(&peak, sizeof(peak), 1, file);
	fout = fwrite(&curvRatio, sizeof(curvRatio), 1, file);
#ifdef ADD_ATTRIBUTE
	fout = fwrite(&log_hessian_eign1, sizeof(log_hessian_eign1), 1, file);
	fout = fwrite(&log_hessian_eign2, sizeof(log_hessian_eign2), 1, file);
	fout = fwrite(&g_hessian_eign1, sizeof(g_hessian_eign1), 1, file);
	fout = fwrite(&g_hessian_eign2, sizeof(g_hessian_eign2), 1, file);
	fout = fwrite(&log_hessian_det, sizeof(log_hessian_det), 1, file);
	fout = fwrite(&log_hessian_trace, sizeof(log_hessian_trace), 1, file);
	fout = fwrite(&g_hessian_det, sizeof(g_hessian_det), 1, file);
	fout = fwrite(&g_hessian_trace, sizeof(g_hessian_trace), 1, file);
	fout = fwrite(&log_dx, sizeof(log_dx), 1, file);
	fout = fwrite(&log_dy, sizeof(log_dy), 1, file);
	fout = fwrite(&log_ds, sizeof(log_ds), 1, file);
	fout = fwrite(&log_dxx, sizeof(log_dxx), 1, file);
	fout = fwrite(&log_dyy, sizeof(log_dyy), 1, file);
	fout = fwrite(&log_dss, sizeof(log_dss), 1, file);
	fout = fwrite(&log_dxy, sizeof(log_dxy), 1, file);
	fout = fwrite(&log_dxs, sizeof(log_dxs), 1, file);
	fout = fwrite(&log_dys, sizeof(log_dys), 1, file);
	fout = fwrite(&g_dx, sizeof(g_dx), 1, file);
	fout = fwrite(&g_dy, sizeof(g_dy), 1, file);
	fout = fwrite(&g_dxx, sizeof(g_dxx), 1, file);
	fout = fwrite(&g_dyy, sizeof(g_dyy), 1, file);
	fout = fwrite(&g_dxy, sizeof(g_dxy), 1, file);
#endif
	fout = fwrite(&descr, sizeof(descr), 1, file);
	fout = fwrite(&pdf, sizeof(pdf), 1, file);
	fout = fwrite(&spatialIndex, sizeof(spatialIndex), 1, file);
	fout = fwrite(&relevance, sizeof(relevance), 1, file);
}

void Feature::fromFile(FILE * file) {
	size_t fout;
	fout = fread(&x, sizeof(x), 1, file);
	fout = fread(&y, sizeof(y), 1, file);
	fout = fread(&scale, sizeof(scale), 1, file);
	fout = fread(&orientation, sizeof(orientation), 1, file);
	fout = fread(&peak, sizeof(peak), 1, file);
	fout = fread(&curvRatio, sizeof(curvRatio), 1, file);
#ifdef ADD_ATTRIBUTE
	fout = fread(&log_hessian_eign1, sizeof(log_hessian_eign1), 1, file);
	fout = fread(&log_hessian_eign2, sizeof(log_hessian_eign2), 1, file);
	fout = fread(&g_hessian_eign1, sizeof(g_hessian_eign1), 1, file);
	fout = fread(&g_hessian_eign2, sizeof(g_hessian_eign2), 1, file);
	fout = fread(&log_hessian_det, sizeof(log_hessian_det), 1, file);
	fout = fread(&log_hessian_trace, sizeof(log_hessian_trace), 1, file);
	fout = fread(&g_hessian_det, sizeof(g_hessian_det), 1, file);
	fout = fread(&g_hessian_trace, sizeof(g_hessian_trace), 1, file);
	fout = fread(&log_dx, sizeof(log_dx), 1, file);
	fout = fread(&log_dy, sizeof(log_dy), 1, file);
	fout = fread(&log_ds, sizeof(log_ds), 1, file);
	fout = fread(&log_dxx, sizeof(log_dxx), 1, file);
	fout = fread(&log_dyy, sizeof(log_dyy), 1, file);
	fout = fread(&log_dss, sizeof(log_dss), 1, file);
	fout = fread(&log_dxy, sizeof(log_dxy), 1, file);
	fout = fread(&log_dxs, sizeof(log_dxs), 1, file);
	fout = fread(&log_dys, sizeof(log_dys), 1, file);
	fout = fread(&g_dx, sizeof(g_dx), 1, file);
	fout = fread(&g_dy, sizeof(g_dy), 1, file);
	fout = fread(&g_dxx, sizeof(g_dxx), 1, file);
	fout = fread(&g_dyy, sizeof(g_dyy), 1, file);
	fout = fread(&g_dxy, sizeof(g_dxy), 1, file);
#endif
	fout = fread(&descr, sizeof(descr), 1, file);
	fout = fread(&pdf, sizeof(pdf), 1, file);
	fout = fread(&spatialIndex, sizeof(spatialIndex), 1, file);
	fout = fread(&relevance, sizeof(relevance), 1, file);
}

void Feature::toStream(BitOutputStream& out) const {
    out.write((unsigned char *)&x, sizeof(x) * 8);
    out.write((unsigned char *)&y, sizeof(y) * 8);
    out.write((unsigned char *)&scale, sizeof(scale) * 8);
    out.write((unsigned char *)&orientation, sizeof(orientation) * 8);
    out.write((unsigned char *)&peak, sizeof(peak) * 8);
    out.write((unsigned char *)&curvRatio, sizeof(curvRatio) * 8);
    out.write((unsigned char *)&descr, sizeof(descr) * 8);
    out.write((unsigned char *)&pdf, sizeof(pdf) * 8);
    out.write((unsigned char *)&spatialIndex, sizeof(spatialIndex) * 8);
    out.write((unsigned char *)&relevance, sizeof(relevance) * 8);
}

void Feature::fromStream(BitInputStream& in) {
    in.read((unsigned char *)&x, sizeof(x) * 8);
    in.read((unsigned char *)&y, sizeof(y) * 8);
    in.read((unsigned char *)&scale, sizeof(scale) * 8);
    in.read((unsigned char *)&orientation, sizeof(orientation) * 8);
    in.read((unsigned char *)&peak, sizeof(peak) * 8);
    in.read((unsigned char *)&curvRatio, sizeof(curvRatio) * 8);
    in.read((unsigned char *)&descr, sizeof(descr) * 8);
    in.read((unsigned char *)&pdf, sizeof(pdf) * 8);
    in.read((unsigned char *)&spatialIndex, sizeof(spatialIndex) * 8);
    in.read((unsigned char *)&relevance, sizeof(relevance) * 8);
}

