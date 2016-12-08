default:stencil

stencil:
	nvcc stencilparallel.cu

clean:
	rm -rf ./a.out

