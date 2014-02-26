SubD Notes
==========

# TODO

- Extend the QuadEdgeMesh each step instead of recalculating it.
- Make original face edges bold in wireframe
- Textures
- Texture coordinate smoothing
- Animate between levels

# Edge based wireframe code

```javascript
this.makeWireGeometry = function() {
    var qe = this.qe;
    var edges = qe.edges;
    var edgeCount = edges.length;
    var verts = this.verts;

    var geometry = new THREE.Geometry2(edgeCount * 4);
    var vertices = geometry.vertices;
    var offset = 0;
    var vert;

    for (var i = 0; i < edgeCount; i++) {
        var edge = edges[i];
        vert = verts[edge.vert0];
        vertices[offset++] = vert.x;
        vertices[offset++] = vert.y;
        vertices[offset++] = vert.z;

        vert = verts[edge.vert1];
        vertices[offset++] = vert.x;
        vertices[offset++] = vert.y;
        vertices[offset++] = vert.z;
    }

    return geometry;
}
```

